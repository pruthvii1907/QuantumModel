use nalgebra_glm as glm;
use glm::Vec3;

#[derive(Clone)]
pub struct Particle{
    pub position: glm::Vec3,
    pub vel: glm::Vec3,
    pub color: [f32; 4],

    }

impl Particle{
    pub fn new(pos: Vec3, color: [f32;4]) -> Self{
        Self{
            position: pos,
            color: color,
            vel: glm::vec3(0.0,0.0,0.0),
            }
        }
    }
