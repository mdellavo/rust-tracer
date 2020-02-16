use std::f64;

use std::sync::Arc;

use crate::utils::Vec3f;
use crate::ray::Ray;
use crate::material::{Material, Hit};

pub trait Hittable: Sync {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<Hit>;
}


#[derive(Clone)]
pub struct Sphere {
    pub center: Vec3f,
    pub radius: f64,
    pub material: Arc<dyn Material>,
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
