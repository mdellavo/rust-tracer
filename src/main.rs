// https://raytracing.github.io/books/RayTracingInOneWeekend.html

extern crate image;
extern crate threadpool;
extern crate indicatif;

use std::f64;
use std::sync::mpsc::channel;
use std::sync::Arc;

use image::{ImageBuffer, Rgba, RgbaImage};
use threadpool::ThreadPool;
use indicatif::{ProgressBar, ProgressStyle};

mod camera;
mod geometry;
mod material;
mod ray;
mod texture;
mod utils;

use camera::Camera;
use geometry::{Hittable, Sphere};
use material::{Dielectric, Diffuse, DiffuseLight, Hit, Material, Metal};
use ray::Ray;
use texture::{CheckerTexture, ConstantTexture, ImageTexture};
use utils::Vec3f;

const WIDTH: u32 = 1024;
const HEIGHT: u32 = 512;

const NUM_SAMPLES: u32 = 250;
const T_MAX: f64 = f64::MAX;

fn color_sky(ray: &Ray) -> Vec3f {
    return Vec3f::new(0.0, 0.0, 0.0);
    let unit = ray.direction().normalize();
    let t = 0.5 * (unit.y + 1.0);
    let c = (1.0 - t) * Vec3f::new(1.0, 1.0, 1.0) + t * Vec3f::new(0.8, 0.8, 1.0);
    return c;
}

fn color(ray: &Ray, world: &impl Hittable, depth: u32) -> Vec3f {
    let result = world.hit(ray, 0.000001, T_MAX);
    match result {
        None => {
            return color_sky(ray);
        }
        Some(hit) => {
            if depth < 50 {
                let result = hit.material.scatter(ray, &hit);
                let emitted = hit.material.emitted(hit.u, hit.v, &hit);
                if result.is_some() {
                    let (attenuation, scattered) = result.unwrap();
                    return emitted
                        + color(&scattered, world, depth + 1).component_mul(&attenuation);
                } else {
                    return emitted;
                }
            }
            return Vec3f::new(0., 0., 0.);
        }
    }
}

fn sample(camera: &Camera, x: u32, y: u32, world: &impl Hittable) -> Vec3f {
    let mut c = Vec3f::new(0., 0., 0.);
    for _ in 0..NUM_SAMPLES {
        let u = ((x as f64) + rand::random::<f64>()) / WIDTH as f64;
        let v = ((y as f64) + rand::random::<f64>()) / HEIGHT as f64;
        let ray = camera.get_ray(u, v);
        c += color(&ray, world, 0);
    }
    c /= NUM_SAMPLES as f64;
    return c;
}

#[derive(Clone)]
struct World<T: Hittable> {
    objects: Vec<T>,
}

impl<T: Hittable> Hittable for World<T> {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<Hit> {
        let mut rv: Option<Hit> = None;
        let mut closest = t_max;
        for hittable in self.objects.iter() {
            let result = hittable.hit(ray, t_min, closest);
            match result {
                Some(hit) => {
                    closest = hit.t;
                    rv = Some(hit);
                }
                None => {}
            }
        }
        return rv;
    }
}

fn random_scenen() -> World<Sphere> {
    let mut objects = Vec::new();

    let ground = Sphere {
        center: Vec3f::new(0., -1000., 0.),
        radius: 1000.,

        material: Arc::new(Diffuse {
            albedo: Arc::new(CheckerTexture {
                odd: Arc::new(ConstantTexture {
                    color: Vec3f::new(0.0, 0.0, 0.0),
                }),
                even: Arc::new(ConstantTexture {
                    color: Vec3f::new(1.0, 1.0, 1.0),
                }),
            }),
        }),
    };
    objects.push(ground);

    let globemap = image::open("earthmap.jpg").unwrap().to_rgba();

    let globe = Sphere {
        center: Vec3f::new(0., 1.0, 0.),
        radius: 0.5,
        material: Arc::new(Diffuse {
            albedo: Arc::new(ImageTexture {
                image: globemap,
                width: WIDTH as f64,
                height: HEIGHT as f64,
            }),
        }),
    };
    objects.push(globe);

    let radius = 5;
    for a in -radius..radius {
        for b in -radius..radius {
                let center = Vec3f::new(a as f64, 0.2, b as f64);
                let sphere;

                match rand::random::<f64>() {
                    0.0..=0.4 => {
                        sphere = Sphere {
                            center: center,
                            radius: 0.2,
                            material: Arc::new(DiffuseLight {
                                texture: Arc::new(ConstantTexture {
                                    color: Vec3f::new_random(),
                                }),
                            }),
                        };
                    }
                    0.4..=0.8 => {
                        sphere = Sphere {
                            center: center,
                            radius: 0.2,
                            material: Arc::new(Metal {
                                albedo: Vec3f::new_random(),
                                fuzz: rand::random::<f64>(),
                            }),
                        };
                    }
                    _ => {
                        sphere = Sphere {
                            center: center,
                            radius: 0.2,
                            material: Arc::new(Dielectric {
                                ref_idx: rand::random::<f64>(),
                            }),
                        };
                    }
                }

                objects.push(sphere);
            }
    }

    return World { objects: objects };
}

fn render_parallel(img: &mut RgbaImage, world: &World<Sphere>, camera: &Camera) {
    let (tx, rx) = channel();
    let pool = ThreadPool::new(8);

    let bar = ProgressBar::new((WIDTH * HEIGHT).into());
    bar.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed} / {eta}] {wide_bar} {percent}%")
            .progress_chars("##-"),
    );

    for x in 0..WIDTH {
        let tx = tx.clone();
        let w = world.clone();
        let c = camera.clone();
        let bar = bar.clone();
        pool.execute(move || {
            for y in 0..HEIGHT {
                tx.send((x, HEIGHT - y - 1, sample(&c, x, y, &w)))
                    .expect("could not dispatch");
                bar.inc(1);
            }
        });
    }
    pool.join();
    drop(tx);

    for (x, y, c) in rx.iter() {
        *img.get_pixel_mut(x, y) = Rgba([
            (255. * c.x.sqrt()) as u8,
            (255. * c.y.sqrt()) as u8,
            (255. * c.z.sqrt()) as u8,
            255,
        ]);
    }
    bar.finish();
}

fn main() {
    let mut img: RgbaImage = ImageBuffer::new(WIDTH, HEIGHT);

    let look_from = Vec3f::new(2.0, 2.0, 2.0);
    let look_at = Vec3f::new(0.0, 0.0, 0.0);
    let up = Vec3f::new(0.0, 1.0, 0.0);
    let fov = 60.0;
    let aspect = WIDTH as f64 / HEIGHT as f64;
    let aperture = 0.1;
    let dist_to_focus = (look_from - look_at).magnitude();
    let camera = Camera::new(look_from, look_at, up, fov, aspect, aperture, dist_to_focus);

    let world = random_scenen();
    render_parallel(&mut img, &world, &camera);

    img.save("out.bmp").unwrap();
    print!("\x07");
}
