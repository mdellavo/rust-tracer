use std::f64;

use crate::utils::Vec3f;

#[derive(Debug)]
pub struct Ray {
    pub a: Vec3f,
    pub b: Vec3f,
}

impl Ray {
    pub fn origin(&self) -> &Vec3f {
        return &self.a;
    }

    pub fn direction(&self) -> &Vec3f {
        return &self.b;
    }

    pub fn point_at(&self, t: f64) -> Vec3f {
        return self.a + (t * self.b);
    }
}
