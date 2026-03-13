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
 
   
fn sample_r(&self, rng: &mut impl Rng) -> f32{
        
        let u: f64 = rng.r#gen();       //creates a random number between 0&1
        
        let idx: usize = cdf_sample(u,&self.cdf_r);
        let r =  idx as f64 * self.r_max / (Self::N_R - 1) as f64;

        r as f32

        }

fn sample_theta(&self, rng: &mut impl Rng) -> f32{

        let u: f64 = rng.r#gen();
        
        let idx: usize = cdf_sample(u,&self.cdf_theta);
        let theta = idx as f64 * PI / (Self::N_THETA - 1) as f64;
        theta as f32
       }


fn sample_phi(&self, rng: &mut impl Rng) -> f32{
//the simplest one
//always uniform - no cdf sample      

        let u: f64 = rng.r#gen();
        
        let res = u * 2.0 * PI;
        res as f32

       }
}


fn probability_flow(pos: &glm::Vec3, m: i32) -> Vec3{
    let r = glm::length(&pos);
    if r < 1e-6 {
        return glm::vec3(0.0, 0.0, 0.0);
        }
    let theta = (pos.y/r).acos();
    let phi = pos.z.atan2(pos.x);
    let sin_theta = theta.sin().max(1e-4); //clamp to 1e-4
    let v_mag = HBAR as i32 * m/ ((M_E * r * sin_theta) as i32);

    glm::vec3((-v_mag * phi.sin() as i32) as f32, 0.0, (v_mag*phi.cos() as i32) as f32)

    }

fn color
/*======================Main============================*/
fn main() {
 
 let pdf = vec![1.0, 2.0, 3.0, 4.0];
    let cdf = build_cdf(&pdf);
    
    for val in &cdf {
        println!("{}", val);
    } 
    println!("{}", spherical_to_cartesian(1.0,0.0,0.0));
    println!("{}", cdf_sample(0.5,&cdf));
}
