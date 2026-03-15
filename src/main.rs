const a0:f64 = 1.0;       //bohr's radius
const HBAR: f32 = 1.0;    //reduced plancks constant
const M_E: f32 = 1.0;     //electron mass

mod octree;
mod particle;

use crate::octree::Octree;
use crate::particle::Particle;

use std::f64::consts::PI;
use rand::thread_rng;
use rand::Rng;
use rand::rngs::ThreadRng;
use nalgebra_glm as glm;
use glm::Vec3;
use glm::vec3;


/*=============BEGINING OF MATH FUNCS=======================*/

/*===========Gamma func==================*/
fn gamma(n:i32)-> f64{
    //gamma for +ve int is n-1 !
    (1..n).map(|k| k as f64).product()
    }


/*===========Laguerre and legendre==================*/
fn assoc_laguerre(k:i32, alpha:i32, x:f64) -> f64{
    if k==0 {return 1.0};
    let mut lag:f64 = 1.0;
    let mut lm1:f64 = 1.0 + alpha as f64 -x;
    if k==1  {lag = lm1;}
    else if k>1 {
    let mut lm2 = 1.0;
    for j in 2..=k
    {
        lag = ((2*j -1 +alpha)as f64 -x) * lm1 - (j-1+alpha) as f64 *lm2;
        lag /= j as f64;
        lm2 = lm1;
        lm1 = lag;
        
    }
   }
   lag
}


//m_abs is the magnetic quantum number 
fn assoc_legendre(l: i32, m_abs: i32, x: f64) -> f64{ 
    let mut pmm: f64 = 1.0;                           
    if m_abs>0 {                                     
        let somx2: f64 = ((1.0-x)*(1.0+x)).sqrt();    
        let mut fact:f64 = 1.0;                       
        for _j in 1..=m_abs{
            pmm *= -fact * somx2;
            fact += 2.0
            }
        }
        let plm: f64 ;
        if l==m_abs {return pmm;}
        else{
            let mut pm1m: f64 = x * ((2 * m_abs + 1) as f64) * pmm;
            if l == m_abs + 1{
                plm = pm1m;
                 }
        else{
            let mut pll: f64 ;
            for ll in m_abs+2..=l{
                pll = ((2 * ll -1) as f64 * x * pm1m - (ll + m_abs - 1) as f64 * pmm) / (ll - m_abs) as f64;
                pmm = pm1m;
                pm1m = pll;

                }
                plm = pm1m;
          }
          }
          plm
    }
/*======================pdf==============================*/
 fn radial_pdf(r:f64, n: i32, l:i32 ) -> f64{
    

        let rho: f64 = 2.0 * r/ ((n as f64) * a0);
        let k: i32 = n - l -1;
        let alpha: i32 = 2 * l + 1;

        let Lag: f64 = assoc_laguerre(k ,alpha,rho);

       let norm: f64 = (2.0 / ((n as f64) * a0)).powi(3) *  gamma(n - l)  /  (2.0 * (n as f64) * gamma(n + l + 1));

       let R: f64 = norm.sqrt()  *  (-rho/2.0).exp()  *  rho.powf(l as f64)  *  Lag; 
       let pdf: f64 = r * r  *  R * R;
        
       pdf 
        
 }
fn angular_pdf(theta:f64, l:i32, m:i32) -> f64{
    let x: f64 = theta.cos();
    let plm: f64 = assoc_legendre(l, m.abs(), x);

    let result:f64 = theta.sin() * plm * plm;

    result

    }
/*===========cdf========================*/
fn build_cdf(pdf_vals: &[f64]) -> Vec<f64>{
    let mut cdf: Vec<f64> = Vec::new();
    let mut sum: f64 = 0.0;
    

//create a cdf number line
    for &vals in pdf_vals {
        sum += vals;
        cdf.push(sum);
        }

//normalize by dividing each value with total
    for v in &mut cdf{
        
        *v /= sum;
        }


    cdf 
    }
/*============convert spherical to cartesian===============*/
fn spherical_to_cartesian(r: f32, theta: f32, phi: f32) -> glm::Vec3{
        let x: f32 = r * theta.sin() * phi.cos();
        let y: f32 = r * theta.cos();
        let z: f32 = r * theta.sin() * phi.sin();
        
        glm::vec3(x,y,z)

        }
/*========random num to physical pos================*/
fn cdf_sample(u:f64, cdf: &[f64]) -> usize {
   
//go to each num & see if num < u yes- go on, no- return idx  

        let idx: usize = cdf.partition_point(|&v| v < u);               
        idx.min(cdf.len() - 1) 

//clamp the result to prevent crash

}
/*=============END OF MATH FUNCS=======================*/

struct quantum_sampler{
    n: i32,
    l: i32,
    m: i32,
    cdf_r: Vec<f64>,
    cdf_theta: Vec<f64>,
    r_max: f64,
    rng: ThreadRng,
    }
 
impl quantum_sampler{

        const N_R: usize = 4096;
        const N_THETA: usize = 2048;

    fn new(n: i32,l: i32,m: i32) -> Self{

        let mut pdf_vals_r: Vec<f64> = Vec::new();
        let mut pdf_vals_a: Vec<f64> = Vec::new();

        let r_max = 10.0 * (n * n) as f64 * a0;
        let dr: f64 = r_max / (Self::N_R-1) as f64;

        for i in 0..Self::N_R{

            let r: f64 = (i as f64) * dr;
            let vals = radial_pdf(r,n,l);
            pdf_vals_r.push(vals);
            
            }
        let cdf_r = build_cdf(&pdf_vals_r);

        let dtheta: f64 = PI / (Self::N_THETA - 1) as f64;

        for j in 0..Self::N_THETA{
                       
            let theta: f64 = j as f64 * dtheta;
            let vals_a = angular_pdf(theta, l, m);
            pdf_vals_a.push(vals_a);

            }
         let cdf_theta = build_cdf(&pdf_vals_a);

        Self{
            n,
            l,
            m,
            r_max,
            cdf_r,
            cdf_theta,
            rng: thread_rng(),
            }
    }
 
   
fn sample_r(&mut self) -> f32{
        
        let u: f64 = self.rng.r#gen();       //creates a random number between 0&1
        
        let idx: usize = cdf_sample(u,&self.cdf_r);
        let r =  idx as f64 * self.r_max / (Self::N_R - 1) as f64;

        r as f32

        }

fn sample_theta(&mut self) -> f32{

        let u: f64 = self.rng.r#gen();
        
        let idx: usize = cdf_sample(u,&self.cdf_theta);
        let theta = idx as f64 * PI / (Self::N_THETA - 1) as f64;
        theta as f32
       }


fn sample_phi(&mut self) -> f32{
//the simplest one
//always uniform - no cdf sample      

        let u: f64 = self.rng.r#gen();
        
        let res = u * 2.0 * PI;
        res as f32

       }

pub fn sample(&mut self) -> Vec3{
    let r = self.sample_r();
    let theta = self.sample_theta();
    let phi = self.sample_phi();

    spherical_to_cartesian(r,theta,phi)
    }
}


fn probability_flow(pos: &glm::Vec3, m: i32) -> Vec3{
    let r = glm::length(&pos);
    if r < 1e-6 {
        return glm::vec3(0.0, 0.0, 0.0);
        }
    let rho = (pos.x * pos.x + pos.y * pos.y).sqrt().max(1e-6);
    let v_mag = HBAR * m as f32 / (M_E * rho as f32);

    let dir = glm::normalize(&glm::vec3(-pos.y, pos.x, 0.0));
    dir * (v_mag as f32)

    }

fn particle_color(psi_squared: f64, max_psi: f64) -> glm::Vec3{
    let t = (psi_squared / max_psi).clamp(0.0, 1.0);
    let mut color: Vec3 = vec3(0.0, 0.0, 0.0);

    if t < 0.25{
        let local_t = t/0.25;
        let black = vec3(0.0,0.0,0.0);
        let purple = vec3(0.5, 0.0, 0.5);
        color = black + (purple - black) * local_t as f32 //purple and black
        }
    else if t < 0.5{
        let local_t = (t-0.25)/0.25;
        let purple = vec3(0.5, 0.0, 0.5);
        let magenta = vec3(1.0, 0.0, 1.0);
        color = purple + (magenta - purple) * local_t as f32;   //magenta and purple
        }
    else if t < 0.75{
        let local_t = (t-0.5)/0.25;
        let magenta = vec3(1.0, 0.0, 1.0);
        let ochre = vec3(0.8, 0.5, 0.0);
        color = magenta + (ochre - magenta) * local_t as f32; //ochre and magenta
        }
    else {
        let local_t = (t-0.75)/0.25;
        let ochre = vec3(0.8, 0.5, 0.0);
        let white = vec3(1.0,1.0,1.0);
        color = ochre + (white - ochre) * local_t as f32; //white
        }

        color
    }

fn particle_gen(n: i32, m:i32, l: i32, num_particles: usize) -> Vec<Particle> {
    
    let mut sampler = quantum_sampler::new(n,l,m);
    let mut particles: Vec<Particle> = Vec::new();

    let mut max_psi: f64 = 0.0;
    for _ in 0..1000{
        let pos = sampler.sample();  //different random positions
        
        //find r, theta and phi
        let r = glm::length(&pos);
        let theta = (pos.z/r).acos();
        let phi = pos.y.atan2(pos.x);

        let psi_sq = radial_pdf(r as f64, n, l) * angular_pdf(theta as f64, l, m);

        if psi_sq > max_psi{
            max_psi = psi_sq;
            }
        }

     for _ in 0..num_particles{
        let pos = sampler.sample(); //another set of random positions
        
        //find r,theta and phi again
        let r = glm::length(&pos);
        let theta = (pos.z/r).acos();
        let phi = pos.y.atan2(pos.x);

        let psi_sq = radial_pdf(r as f64,n,l) * angular_pdf(theta as f64,l, m);
        let color = particle_color(psi_sq,max_psi);
        let velocity = probability_flow(&pos, m);

        let color_array = [color.x, color.y, color.z, 1.0];
        let particle = Particle::new(pos,color_array);

        particles.push(particle);
         }

    particles
    }
/*======================Camera============================*/
struct Camera {
    position: Vec3,
    target: Vec3,
    up: Vec3,
    fov: f32,   //field of view
    aspect: f32,
    near: f32,
    far: f32,

    //Orbital camera params
    distance: f32,  //distance from target
    theta: f32,     //vertical angle
    phi: f32,       //horizontal angle
    }

impl Camera{
    pub fn new(position: Vec3, target: Vec3, up: Vec3, aspect: f32) -> Self{
        let to_target = target - position;
        let distance = glm::length(&to_target);

        let theta = (to_target.y/distance).acos();
        let phi = to_target.z.atan2(to_target.x);

        Self{
            position,
            target,
            up,
            fov: 45.0,
            aspect,
            near: 0.1,
            far: 1000.0,
            distance,
            theta,
            phi,
            }
}

//view_matrix to tell the camera where it is looking
pub fn view_matrix(&self) -> glm::Mat4{
    glm::look_at(&self.position, &self.up, &self.target)
    }

//projection_matrix to gain perspective make far away things smaller
pub fn projection_matrix(&self) -> glm::Mat4{
    glm::perspective(self.aspect, self.fov.to_radians(), self.near, self.far)
    }

pub fn orbit(&mut self, delta_theta: f32, delta_phi: f32){ //change the theta and phi calues and then calculate new pos
    self.theta += delta_theta;
    self.phi += delta_phi;

    self.theta = self.theta.clamp(0.1, std::f32::consts::PI - 0.1);  //clamp the theta to 0.1 to pi-1 to prevent flipping

    let x = self.distance * self.theta.sin() * self.phi.cos();   //calculate new x,y and z values
    let y = self.distance * self.theta.cos();
    let z = self.distance * self.theta.sin() * self.phi.sin();

    self.position = self.target + glm::vec3(x,y,z);         //update the position
 }

pub fn zoom(&mut self, delta: f32) {
    //same as orbit find clamp update
    self.distance += delta;
    
    self.distance = self.distance.clamp(1.0, 500.0);
    
    let x = self.distance * self.theta.sin() * self.phi.cos();
    let y = self.distance * self.theta.cos();
    let z = self.distance * self.theta.sin() * self.phi.sin();
    }
}
/*===========================================Main==================================================*/
fn main() {
    println!("Generating particles for 2p orbital (n=2, l=1, m=0)...");
    
    let particles = particle_gen(2, 0, 1, 1000);  // n=2, m=0, l=1, 1000 particles
    
    println!("Generated {} particles", particles.len());
    
    // Print first 5 particles to see what they look like
    for (i, particle) in particles.iter().take(5).enumerate() {
        println!("Particle {}: pos = ({:.3}, {:.3}, {:.3}), color = ({:.3}, {:.3}, {:.3}, {:.3})", 
            i,
            particle.position.x, particle.position.y, particle.position.z,
            particle.color[0], particle.color[1], particle.color[2], particle.color[3]
        );
    }
} 
