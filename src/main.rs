// https://raytracing.github.io/books/RayTracingInOneWeekend.html

extern crate image;
extern crate nalgebra;
extern crate threadpool;

use std::sync::mpsc::channel;

use threadpool::ThreadPool;
use image::{ImageBuffer, Rgba, RgbaImage};
use nalgebra::Vector3;

const WIDTH :u32 = 1024;
const HEIGHT :u32 = 512;

const NUM_SAMPLES :u32 = 25;


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


trait Hittable {
    fn hit(&self, ray: &Ray) -> f64;
}

#[derive(Copy, Clone)]
struct Sphere {
  center: Vec3f,
  radius: f64,
}


impl Hittable for Sphere {
    fn hit(&self, ray: &Ray) -> f64 {
        let oc = ray.origin() - self.center;
        let a = ray.direction().dot(ray.direction());
        let b = 2. * oc.dot(ray.direction());
        let c = oc.dot(&oc) - (self.radius * self.radius);
        let discriminant = (b * b) - (4. * a * c);
        if discriminant < 0.0 {
            return -1.0;
        } else {
            return (-b - discriminant.sqrt()) / (2.0 * a);
        }
    }
}

fn color_sky(ray: &Ray) -> Vec3f {
    let unit = ray.direction().normalize();
    let t = 0.5 * (unit.y + 1.0);
    let c = (1.0 - t) * Vec3f::new(1.0, 1.0, 1.0) + t * Vec3f::new(0.5, 0.7, 1.0);
    return c;
}

fn color(ray: &Ray) -> Vec3f {
    let sphere = Sphere {
        center: Vec3f::new(0., 0., -1.),
        radius: 0.5,
    };
    let t = sphere.hit(ray);
    if  t > 0.0 {
        let n = (ray.point_at(t) - Vec3f::new(0.0, 0.0, -1.0)).normalize();
        return 0.5 * Vec3f::new(n.x + 1., n.y + 1., n.z + 1.);
    }

    return color_sky(ray);
}

fn sample(camera: &Camera, x: u32, y: u32) -> Vec3f {
    let mut c = Vec3f::new(0., 0., 0.);
    for _ in 0..NUM_SAMPLES {
        let u = ((x as f64) + rand::random::<f64>()) / WIDTH as f64;
        let v = ((y as f64) + rand::random::<f64>()) / HEIGHT as f64;
        let ray = camera.get_ray(u, v);
        c += color(&ray);
    }
    c /= NUM_SAMPLES as f64;
    return c;
}


fn main() {
    let mut img : RgbaImage = ImageBuffer::new(WIDTH, HEIGHT);

    let camera = Camera {
        lower_left_corner: Vec3f::new(-2.0, -1.0, -1.0),
        horizontal: Vec3f::new(4.0, 0.0, 0.0),
        vertical: Vec3f::new(0.0, 2.0, 0.0),
        origin: Vec3f::new(0.0, 0.0, 0.0),
    };

    let (tx, rx) = channel();
    let pool = ThreadPool::new(8);
    for x in 0..WIDTH {
        for y in 0..HEIGHT {
            let tx = tx.clone();
            pool.execute(move || {
                tx.send((x, HEIGHT - y - 1, sample(&camera, x, y)))
                    .expect("could not dispatch");
            });
        }
    }
    drop(tx);
    pool.join();

    for (x, y, c) in rx.iter() {
        *img.get_pixel_mut(x, y) = Rgba([
            (255. * c.x) as u8,
            (255. * c.y) as u8,
            (255. * c.z) as u8,
            255,
        ]);
    }

    img.save("out.png").unwrap();
}
