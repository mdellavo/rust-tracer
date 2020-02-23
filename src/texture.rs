use std::f64;

use std::sync::Arc;

use image::RgbaImage;

use crate::utils::Vec3f;
use crate::material::Hit;

pub trait Texture: Sync + Send {
    fn value(&self, u: f64, v: f64, p: &Hit) -> Vec3f;
}

pub struct ConstantTexture {
    pub color: Vec3f,
}

impl Texture for ConstantTexture {
    fn value(&self, u: f64, v: f64, p: &Hit) -> Vec3f {
        return self.color;
    }
}


pub struct CheckerTexture {
    pub odd: Arc<dyn Texture>,
    pub even: Arc<dyn Texture>,
}

impl Texture for CheckerTexture {
  fn value(&self, u: f64, v: f64, p: &Hit) -> Vec3f {
        let sines = (10.0 * p.p.x).sin() * (10.0 * p.p.y).sin() * (10.0 * p.p.z).sin();
        if sines < 0.0 {
            return self.odd.value(u, v, p);
        } else {
            return self.even.value(u, v, p);
        }
    }
}


pub struct ImageTexture {
    pub image: RgbaImage,
    pub width: f64,
    pub height: f64,
}

impl Texture for ImageTexture {
    fn value(&self, u: f64, v: f64, p: &Hit) -> Vec3f {
        let mut i = u * self.width;
        let mut j = (1.0-v) * self.height - 0.001;

        if i < 0.0 {
            i = 0.0;
        }
        if i > self.width - 1.0 {
            i = self.width -1.0;
        }

        if j < 0.0 {
            j = 0.0;
        }
        if j > self.height - 1.0 {
            j = self.height -1.0;
        }

        let pixel = self.image.get_pixel(i as u32, j as u32);
        let rgba = *pixel;
        return Vec3f::new(rgba[0] as f64 / 255.0, rgba[1] as f64 / 255.0, rgba[2] as f64 / 255.0);
    }
}
