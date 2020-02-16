use std::f64;

use crate::utils::{random_point_in_unit, Vec3f};
use crate::ray::Ray;


#[derive(Copy, Clone)]
pub struct Camera {
    pub lower_left_corner: Vec3f,
    pub horizontal: Vec3f,
    pub vertical: Vec3f,
    pub origin: Vec3f,
    pub u: Vec3f,
    pub v: Vec3f,
    pub w: Vec3f,
    pub lense_radius: f64,
}

impl Camera {
    pub fn new(
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

    pub fn get_ray(&self, s: f64, t: f64) -> Ray {
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
