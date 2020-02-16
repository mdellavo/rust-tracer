// https://raytracing.github.io/books/RayTracingInOneWeekend.html

extern crate image;
extern crate nalgebra;
extern crate threadpool;

use std::f64;
use std::sync::mpsc::channel;
use std::sync::Arc;

use image::{ImageBuffer, Rgba, RgbaImage};
use nalgebra::Vector3;
use threadpool::ThreadPool;

const WIDTH: u32 = 1024;
const HEIGHT: u32 = 512;

const NUM_SAMPLES: u32 = 250;
const T_MAX: f64 = f64::MAX;

type Vec3f = Vector3<f64>;

#[derive(Copy, Clone)]
struct Camera {
    lower_left_corner: Vec3f,
    horizontal: Vec3f,
    vertical: Vec3f,
    origin: Vec3f,
    u: Vec3f,
    v: Vec3f,
    w: Vec3f,
    lense_radius: f64,
}

impl Camera {
    fn new(
        look_from: Vec3f,
        look_at: Vec3f,
        up: Vec3f,
        fov: f64,
        aspect: f64,
        aperture: f64,
        focus_dist: f64,
    ) -> Camera {
        let (u, v, w);

        let theta = fov * f64::consts::PI / 180.0;
        let half_height = (theta / 2.0).tan();
        let half_width = aspect * half_height;
        let lense_radius = aperture / 2.0;

        w = (look_from - look_at).normalize();
        u = up.cross(&w).normalize();
        v = w.cross(&u);

        let origin = look_from;
        let lower_left_corner = origin
            - (half_width * u * focus_dist)
            - (half_height * v * focus_dist)
            - (w * focus_dist);
        let horizontal = 2.0 * half_width * u * focus_dist;
        let vertical = 2.0 * half_height * v * focus_dist;

        return Camera {
            origin: origin,
            lower_left_corner: lower_left_corner,
            vertical: vertical,
            horizontal: horizontal,
            u: u,
            v: v,
            w: w,
            lense_radius: lense_radius,
        };
    }

    fn get_ray(&self, s: f64, t: f64) -> Ray {
        let rd = self.lense_radius * random_point_in_unit();
        let offset = self.u * rd.x + self.v * rd.y;
        let ray = Ray {
            a: self.origin + offset,
            b: self.lower_left_corner + (s * self.horizontal) + (t * self.vertical)
                - self.origin
                - offset,
        };
        return ray;
    }
}

#[derive(Debug)]
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

trait Material: Sync + Send {
    fn scatter(&self, ray: &Ray, hit: &Hit) -> Option<(Vec3f, Ray)>;
}

#[derive(Copy, Clone)]
struct Diffuse {
    albedo: Vec3f,
}

impl Material for Diffuse {
    fn scatter(&self, ray: &Ray, hit: &Hit) -> Option<(Vec3f, Ray)> {
        let target = hit.p + hit.normal + random_point_in_unit();
        let scattered = Ray {
            a: hit.p,
            b: target - hit.p,
        };
        return Some((self.albedo.clone(), scattered));
    }
}

fn reflect(a: &Vec3f, b: &Vec3f) -> Vec3f {
    return a - 2.0 * a.dot(b) * b;
}

struct Metal {
    albedo: Vec3f,
    fuzz: f64,
}

impl Material for Metal {
    fn scatter(&self, ray: &Ray, hit: &Hit) -> Option<(Vec3f, Ray)> {
        let reflected = reflect(&ray.direction().normalize(), &hit.normal);
        let scattered = Ray {
            a: hit.p,
            b: reflected + (self.fuzz * random_point_in_unit()),
        };

        if scattered.direction().dot(&hit.normal) > 0. {
            return Some((self.albedo.clone(), scattered));
        }

        return None;
    }
}

fn refract(v: &Vec3f, n: &Vec3f, ni_over_nt: f64) -> Option<Vec3f> {
    let uv = v.normalize();
    let dt = uv.dot(&n);
    let discriminant = 1.0 - ni_over_nt * ni_over_nt * (1.0 - dt * dt);
    if discriminant > 0.0 {
        let refracted = ni_over_nt * (uv - n * dt) - n * discriminant.sqrt();
        return Some(refracted);
    }
    return None;
}

fn schlick(cosine: f64, ref_idx: f64) -> f64 {
    let mut r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
    r0 = r0 * r0;
    return r0 + (1.0 - r0) * (1.0 - cosine).powi(5);
}

struct Dielectric {
    ref_idx: f64,
}

impl Material for Dielectric {
    fn scatter(&self, ray: &Ray, hit: &Hit) -> Option<(Vec3f, Ray)> {
        let reflected = reflect(&ray.direction(), &hit.normal);
        let attenuation = Vec3f::new(1.0, 1.0, 1.0);

        let outward_normal;
        let ni_over_nt;
        let cosine: f64;
        let reflect_prob: f64;
        if ray.direction().dot(&hit.normal) > 0.0 {
            outward_normal = -hit.normal.clone();
            ni_over_nt = self.ref_idx;
            cosine = self.ref_idx * ray.direction().dot(&hit.normal) / ray.direction().magnitude();
        } else {
            outward_normal = hit.normal.clone();
            ni_over_nt = 1.0 / self.ref_idx;
            cosine = -ray.direction().dot(&hit.normal) / ray.direction().magnitude();
        }

        let result = refract(&ray.direction(), &outward_normal, ni_over_nt);
        if result.is_some() {
            reflect_prob = schlick(cosine, self.ref_idx);
        } else {
            reflect_prob = 1.0;
        }

        let scattered;
        if rand::random::<f64>() < reflect_prob {
            scattered = Ray {
                a: hit.p,
                b: reflected,
            };
        } else {
            scattered = Ray {
                a: hit.p,
                b: result.unwrap(),
            };
        }
        return Some((attenuation, scattered));
    }
}

#[derive(Clone)]
struct Sphere {
    center: Vec3f,
    radius: f64,
    material: Arc<dyn Material>,
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
                    material: self.material.clone(),
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

struct Hit {
    t: f64,
    p: Vec3f,
    normal: Vec3f,
    material: Arc<dyn Material>,
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
            albedo: Vec3f::new(0.5, 0.5, 0.5),
        }),
    };
    objects.push(ground);

    for a in -11..11 {
        for b in -11..11 {
            let center = Vec3f::new(
                2.0 * a as f64 * rand::random::<f64>(),
                0.2,
                2.0 * b as f64 * rand::random::<f64>(),
            );
            let sphere;

            match rand::random::<f64>() {
                0.0..=0.4 => {
                    sphere = Sphere {
                        center: center,
                        radius: 0.2,
                        material: Arc::new(Diffuse {
                            albedo: Vec3f::new_random(),
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
    let aperture = 0.0;
    let dist_to_focus = (look_from - look_at).magnitude();
    let camera = Camera::new(look_from, look_at, up, fov, aspect, aperture, dist_to_focus);

    let world = random_scenen();
    render_parallel(&mut img, &world, &camera);

    img.save("out.bmp").unwrap();
}
