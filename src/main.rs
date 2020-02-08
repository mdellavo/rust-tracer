// https://raytracing.github.io/books/RayTracingInOneWeekend.html

extern crate image;
extern crate nalgebra;
extern crate threadpool;

use std::f64;
use std::sync::mpsc::channel;

use image::{ImageBuffer, Rgba, RgbaImage};
use nalgebra::Vector3;
use threadpool::ThreadPool;

const WIDTH: u32 = 1024;
const HEIGHT: u32 = 512;

const NUM_SAMPLES: u32 = 5;
const T_MAX :f64 = f64::MAX;

type Vec3f = Vector3<f64>;

#[derive(Copy, Clone)]
struct Camera {
    lower_left_corner: Vec3f,
    horizontal: Vec3f,
    vertical: Vec3f,
    origin: Vec3f,
}

impl Camera {
    fn get_ray(&self, u: f64, v: f64) -> Ray {
        let ray = Ray {
            a: self.origin,
            b: self.lower_left_corner + (u * self.horizontal) + (v * self.vertical),
        };
        return ray;
    }
}

#[derive(Copy, Clone)]
struct Ray {
    a: Vec3f,
    b: Vec3f,
}

impl Ray {
    fn origin(&self) -> &Vec3f {
        return &self.a;
    }

    fn direction(&self) -> &Vec3f {
        return &self.b;
    }

    fn point_at(&self, t: f64) -> Vec3f {
        return self.a + (t * self.b);
    }
}

trait Hittable: Sync {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<Hit>;
}

#[derive(Copy, Clone)]
struct Sphere {
    center: Vec3f,
    radius: f64,
}

impl Hittable for Sphere {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<Hit> {
        let oc = ray.origin() - self.center;
        let a = ray.direction().dot(ray.direction());
        let b = 2. * oc.dot(ray.direction());
        let c = oc.dot(&oc) - (self.radius * self.radius);
        let discriminant = (b * b) - (4. * a * c);
        if discriminant < 0.0 {
            return None;
        }

        let _check = move |t: f64| -> Option<Hit> {
            if t > t_min && t < t_max {
                let p = ray.point_at(t);
                return Some(Hit {
                    t: t,
                    p: p,
                    normal: ((p - self.center) / self.radius).normalize(),
                });
            }
            return None;
        };

        let mut result = _check((-b - discriminant.sqrt()) / (2. * a));
        if result.is_none() {
            result = _check((-b + discriminant.sqrt()) / (2. * a));
        }
        return result;
    }
}

fn color_sky(ray: &Ray) -> Vec3f {
    let unit = ray.direction().normalize();
    let t = 0.5 * (unit.y + 1.0);
    let c = (1.0 - t) * Vec3f::new(1.0, 1.0, 1.0) + t * Vec3f::new(0.5, 0.7, 1.0);
    return c;
}


fn random_point_in_unit() -> Vec3f {
    return Vec3f::new_random().normalize();
}


fn color(ray: &Ray, world: &dyn Hittable) -> Vec3f {
    let result = world.hit(ray, 0.001, T_MAX);
    match result {
        None => {
            return color_sky(ray);
        }
        Some(hit) => {
            let target = hit.p + hit.normal + random_point_in_unit();
            let r = Ray {
                a: hit.p,
                b: target - hit.p
            };
            return 0.5 * color(&r, world);
        }
    }
}

fn sample(camera: &Camera, x: u32, y: u32, world: &dyn Hittable) -> Vec3f {
    let mut c = Vec3f::new(0., 0., 0.);
    for _ in 0..NUM_SAMPLES {
        let u = ((x as f64) + rand::random::<f64>()) / WIDTH as f64;
        let v = ((y as f64) + rand::random::<f64>()) / HEIGHT as f64;
        let ray = camera.get_ray(u, v);
        c += color(&ray, world);
    }
    c /= NUM_SAMPLES as f64;
    return c;
}


struct Hit {
    t: f64,
    p: Vec3f,
    normal: Vec3f,
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
                },
                None => {},
            }
        }
        return rv;
    }
}

fn main() {
    let mut img: RgbaImage = ImageBuffer::new(WIDTH, HEIGHT);

    let camera = Camera {
        lower_left_corner: Vec3f::new(-2.0, -1.0, -1.0),
        horizontal: Vec3f::new(4.0, 0.0, 0.0),
        vertical: Vec3f::new(0.0, 2.0, 0.0),
        origin: Vec3f::new(0.0, 0.0, 0.0),
    };
    let sphere = Sphere {
        center: Vec3f::new(0., 0., -1.),
        radius: 0.5,
    };
    let ground = Sphere {
        center: Vec3f::new(0., -100.5, -1.),
        radius: 100.,
    };

    let world = World {
        objects: vec![sphere, ground],
    };

    let (tx, rx) = channel();
    let pool = ThreadPool::new(8);
    for x in 0..WIDTH {
        for y in 0..HEIGHT {
            let tx = tx.clone();
            let w = world.clone();
            pool.execute(move || {
                tx.send((x, HEIGHT - y - 1, sample(&camera, x, y, &w)))
                    .expect("could not dispatch");
            });
        }
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

    img.save("out.png").unwrap();
}
