use std::f64;

use std::sync::Arc;

use crate::utils::Vec3f;
use crate::material::Hit;

pub trait Texture: Sync + Send {
    fn value(&self, u: f64, v: f64, p: &Hit) -> &Vec3f;
}


pub struct ConstantTexture {
    pub color: Vec3f,
}

impl Texture for ConstantTexture {
    fn value(&self, u: f64, v: f64, p: &Hit) -> &Vec3f {
        return &self.color;
    }
}


pub struct CheckerTexture {
    pub odd: Arc<dyn Texture>,
    pub even: Arc<dyn Texture>,
}

impl Texture for CheckerTexture {
  fn value(&self, u: f64, v: f64, p: &Hit) -> &Vec3f {
        let sines = (10.0 * p.p.x).sin() * (10.0 * p.p.y).sin() * (10.0 * p.p.z).sin();
        if sines < 0.0 {
            return self.odd.value(u, v, p);
        } else {
            return self.even.value(u, v, p);
        }
    }
}
