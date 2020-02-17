use std::f64;
use std::sync::Arc;

use crate::utils::{Vec3f, random_point_in_unit};
use crate::ray::Ray;
use crate::texture::Texture;


pub struct Hit {
    pub t: f64,
    pub p: Vec3f,
    pub normal: Vec3f,
    pub material: Arc<dyn Material>,
}


pub trait Material: Sync + Send {
    fn scatter(&self, ray: &Ray, hit: &Hit) -> Option<(Vec3f, Ray)>;
}

pub struct Diffuse {
    pub albedo: Arc<dyn Texture>,
}

impl Material for Diffuse {
    fn scatter(&self, ray: &Ray, hit: &Hit) -> Option<(Vec3f, Ray)> {
        let target = hit.p + hit.normal + random_point_in_unit();
        let scattered = Ray {
            a: hit.p,
            b: target - hit.p,
        };
        let attenuation = self.albedo.value(0.0, 0.0, &hit);
        return Some((*attenuation, scattered));
    }
}

fn reflect(a: &Vec3f, b: &Vec3f) -> Vec3f {
    return a - 2.0 * a.dot(b) * b;
}

pub struct Metal {
    pub albedo: Vec3f,
    pub fuzz: f64,
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

pub struct Dielectric {
    pub ref_idx: f64,
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
