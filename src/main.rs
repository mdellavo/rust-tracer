// https://raytracing.github.io/books/RayTracingInOneWeekend.html

extern crate image;
extern crate threadpool;

use std::f64;
use std::sync::mpsc::channel;
use std::sync::Arc;

use image::{ImageBuffer, Rgba, RgbaImage};
use threadpool::ThreadPool;

mod utils;
mod camera;
mod material;
mod ray;
mod geometry;
mod texture;

use utils::Vec3f;
use camera::Camera;
use material::{Hit, Material, Dielectric, Diffuse, Metal};
use ray::Ray;
use geometry::{Hittable, Sphere};
use texture::{ConstantTexture, CheckerTexture};

const WIDTH: u32 = 1024;
const HEIGHT: u32 = 512;

const NUM_SAMPLES: u32 = 250;
const T_MAX: f64 = f64::MAX;



fn color_sky(ray: &Ray) -> Vec3f {
    let unit = ray.direction().normalize();
    let t = 0.5 * (unit.y + 1.0);
    let c = (1.0 - t) * Vec3f::new(1.0, 1.0, 1.0) + t * Vec3f::new(0.5, 0.7, 1.0);
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
                if result.is_some() {
                    let (attenuation, scattered) = result.unwrap();
                    return color(&scattered, world, depth + 1).component_mul(&attenuation);
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
                odd:  Arc::new(ConstantTexture { color: Vec3f::new(0.0, 0.0, 0.0) }),
                even:  Arc::new(ConstantTexture { color: Vec3f::new(1.0, 1.0, 1.0)}),
            }),
        }),
    };
    objects.push(ground);

    for a in -25..25 {
        for b in -25..25 {
            let center = Vec3f::new(
                a as f64,//2.0 * a as f64 * rand::random::<f64>(),
                0.2,
                b as f64,//2.0 * b as f64 * rand::random::<f64>(),
            );
            let sphere;

            match rand::random::<f64>() {
                0.0..=0.4 => {
                    sphere = Sphere {
                        center: center,
                        radius: 0.2,
                        material: Arc::new(Diffuse {
                            albedo: Arc::new(ConstantTexture {
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
    for x in 0..WIDTH {
        let tx = tx.clone();
        let w = world.clone();
        let c = camera.clone();
        pool.execute(move || {
            for y in 0..HEIGHT {
                tx.send((x, HEIGHT - y - 1, sample(&c, x, y, &w)))
                    .expect("could not dispatch");
            }
            println!("{} done!", x);
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
}

fn main() {
    let mut img: RgbaImage = ImageBuffer::new(WIDTH, HEIGHT);

    let look_from = Vec3f::new(1.0, 1.0, 1.0);
    let look_at = Vec3f::new(0.0, 0.0, 0.0);
    let up = Vec3f::new(0.0, 1.0, 0.0);
    let fov = 90.0;
    let aspect = WIDTH as f64 / HEIGHT as f64;
    let aperture = 0.01;
    let dist_to_focus = (look_from - look_at).magnitude();
    let camera = Camera::new(look_from, look_at, up, fov, aspect, aperture, dist_to_focus);

    let world = random_scenen();
    render_parallel(&mut img, &world, &camera);

    img.save("out.bmp").unwrap();
}
