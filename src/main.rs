use std::f64::consts::PI;
use num::complex::Complex;
use std::time::Instant;

fn complex_modulo(z1: Complex<f64>, z2: Complex<f64>) -> Complex<f64> {
    // Calculate the absolute values of the complex numbers
    let abs_z1 = z1.norm();
    let abs_z2 = z2.norm();
    
    // Calculate the modulus
    let remainder = abs_z1 % abs_z2;
    
    // Return the remainder as a complex number with the same argument as z1
    let angle_z1 = z1.arg();
    Complex::from_polar(&remainder, &angle_z1)
}

fn sin_complex(v: Complex<f64>) -> Complex<f64> {
    // Take advantage of periodicity of trigonometric functions to evaluate at large values
    let a = v.re % (2.0 * PI);
    let b = complex_modulo(v.im * Complex::i(), 2.0 * PI * Complex::i());

    Complex::new(a.sin() * b.im.cosh(), a.cos() * b.im.sinh())
}

fn cos_complex(v: Complex<f64>) -> Complex<f64> {
    // Take advantage of periodicity of trigonometric functions to evaluate at large values
    let a = v.re % (2.0 * PI);
    let b = complex_modulo(v.im * Complex::i(), 2.0 * PI * Complex::i());

    Complex::new(a.cos() * b.im.cosh(), -a.sin() * b.im.sinh())
}

fn transfer_matrix(result: &mut [[Complex<f64>; 2]; 2], k_0: f64, n: Complex<f64>, d: f64, theta: f64) {
    let k_z = k_0 * n * theta.cos(); // Calculate longitudinal K

    // Reduce redundant calculations in construction T_i
    let q_1 = cos_complex(k_z * d);
    let q_2 = Complex::i() * sin_complex(k_z * d);
    let n_cos_th = n * theta.cos();

    // Transfer matrix calculation derived from Maxwell's equations
    result[0][0] = q_1;
    result[0][1] = q_2 / n_cos_th;
    result[1][0] = n_cos_th * q_2;
    result[1][1] = q_1;
}

fn solve_tmm(
    r: &mut f64,
    t: &mut f64,
    layers: &Vec<[Complex<f64>; 2]>,
    num_layers: usize,
    wavelength: f64,
    theta: f64,
) {
    let k_0 = (2.0 * PI) / wavelength;

    // Global transfer matrix initialization
    let mut m = [[Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
                [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)]];

    let n_0 = layers[0][0];
    let n_l = layers[num_layers - 1][0];

    for i in 1..num_layers - 1 {
        let n = layers[i][0];
        let d = layers[i][1].re;

        let theta_i = (n_0.re * theta.sin() / n.re).asin();
        let mut t_i = [[Complex::default(); 2]; 2];
        transfer_matrix(&mut t_i, k_0, n, d, theta_i);

        // Update global transfer matrix with dot product
        let mut temp = [[Complex::default(); 2]; 2];
        for j in 0..2 {
            for k in 0..2 {
                temp[j][k] = (0..2).map(|l| m[j][l] * t_i[l][k]).sum();
            }
        }
        m = temp;
    }

    let q_1 = n_0 * theta.cos();
    let theta_l = (n_0.re * theta.sin() / n_l.re).asin();
    let q_2 = n_l * theta_l.cos();

    let m_in = [
        [Complex::new(1.0, 0.0), Complex::new(1.0, 0.0)],
        [q_1, -q_1],
    ];

    let m_out = [
        [Complex::new(1.0, 0.0), Complex::new(1.0, 0.0)],
        [q_2, -q_2],
    ];

    let det = m_in[0][0] * m_in[1][1] - m_in[0][1] * m_in[1][0];
    let m_inv = [
        [m_in[1][1] / det, -m_in[0][1] / det],
        [-m_in[1][0] / det, m_in[0][0] / det],
    ];

    let mut m_temp = [[Complex::default(); 2]; 2];
    for j in 0..2 {
        for k in 0..2 {
            m_temp[j][k] = (0..2).map(|l| m[j][l] * m_out[l][k]).sum();
        }
    }

    let mut m_total = [[Complex::default(); 2]; 2];
    for j in 0..2 {
        for k in 0..2 {
            m_total[j][k] = (0..2).map(|l| m_inv[j][l] * m_temp[l][k]).sum();
        }
    }

    let r_complex = m_total[1][0] / m_total[0][0];
    let t_complex = Complex::new(1.0, 0.0) / m_total[0][0];
    *r = r_complex.norm_sqr();
    *t = t_complex.norm_sqr() * (n_l * theta_l.cos()).re / (n_0 * theta.cos()).re;
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() % 2 != 1 {
        eprintln!("Usage: {} <complex double> <integer> ...", args[0]);
        std::process::exit(1);
    }

    let start = Instant::now();

    let mut r = 0.0;
    let mut t = 0.0;
    let num_layers = (args.len() - 3) / 2;
    let mut layers = Vec::new();

    for i in 0..num_layers {
        let real_part: f64 = args[2 * i + 3].parse().unwrap();
        let imag_part: f64 = args[2 * i + 4].parse().unwrap();
        layers.push([Complex::new(real_part, imag_part), args[2 * i + 5].parse().unwrap()]);
    }

    let theta: f64 = args[2].parse().unwrap();
    let wl: f64 = args[1].parse().unwrap();
    solve_tmm(&mut r, &mut t, &layers, num_layers, wl, theta);

    println!("Reflectance: {}", r);
    println!("Transmittance: {}", t);

    let duration = start.elapsed();
    println!("Execution Time: {:.2?} ms", duration.as_micros());
}
