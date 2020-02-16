extern crate nalgebra;

use std::f64;

use nalgebra::Vector3;

pub type Vec3f = Vector3<f64>;

pub fn random_point_in_unit() -> Vec3f {
    return Vec3f::new_random().normalize() - Vec3f::new(1.0, 1.0, 1.0);
}
