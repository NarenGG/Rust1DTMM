use num::complex::Complex;
use std::f64::consts::PI;
use std::time::Instant;

fn complex_modulo(z1: Complex<f64>, z2: Complex<f64>) -> Complex<f64> {
    let abs_z1 = z1.norm();
    let abs_z2 = z2.norm();
    let remainder = abs_z1 % abs_z2;
    let angle_z1 = z1.arg();
    Complex::from_polar(&remainder, &angle_z1)
}

fn sin_complex(v: Complex<f64>) -> Complex<f64> {
    let a = v.re % (2.0 * PI);
    let b = complex_modulo(Complex::new(0.0, v.im), Complex::new(0.0, 2.0 * PI));
    Complex::new(a.sin() * b.im.cosh(), a.cos() * b.im.sinh())
}

fn cos_complex(v: Complex<f64>) -> Complex<f64> {
    let a = v.re % (2.0 * PI);
    let b = complex_modulo(Complex::new(0.0, v.im), Complex::new(0.0, 2.0 * PI));
    Complex::new(a.cos() * b.im.cosh(), -a.sin() * b.im.sinh())
}

fn transfer_matrix(
    result: &mut [[Complex<f64>; 2]; 2],
    k_0: f64,
    n: Complex<f64>,
    d: f64,
    theta: f64,
) {
    let k_z = k_0 * n * theta.cos();
    let q_1 = cos_complex(k_z * d);
    let q_2 = Complex::i() * sin_complex(k_z * d);
    let n_cos_th = n * theta.cos();
    result[0][0] = q_1;
    result[0][1] = q_2 / n_cos_th;
    result[1][0] = n_cos_th * q_2;
    result[1][1] = q_1;
}

fn solve_tmm(
    r: &mut f64,
    t: &mut f64,
    layers: &Vec<[Complex<f64>; 2]>,
    wavelength: f64,
    theta: f64,
) {
    let k_0 = (2.0 * PI) / wavelength;
    let mut m = [[Complex::new(1.0, 0.0); 2]; 2];
    let n_0 = layers[0][0];
    let n_l = layers[layers.len() - 1][0];

    for i in 1..layers.len() - 1 {
        let n = layers[i][0];
        let d = layers[i][1].re;
        let theta_i = (n_0.re * theta.sin() / n.re).asin();
        let mut t_i = [[Complex::<f64>::new(0.0, 0.0); 2]; 2];
        transfer_matrix(&mut t_i, k_0, n, d, theta_i);

        let mut temp = [[Complex::new(0.0, 0.0); 2]; 2];
        for j in 0..2 {
            for k in 0..2 {
                temp[j][k] = Complex::new(0.0, 0.0);
                for l in 0..2 {
                    temp[j][k] += m[j][l] * t_i[l][k];
                }
            }
        }
        m = temp;
    }

    let q_1 = n_0 * theta.cos();
    let theta_l = (n_0.re * theta.sin() / n_l.re).asin();
    let q_2 = n_l * theta_l.cos();
    let m_in = [[Complex::new(1.0, 0.0), Complex::new(1.0, 0.0)], [q_1, -q_1]];
    let m_out = [[Complex::new(1.0, 0.0), Complex::new(1.0, 0.0)], [q_2, -q_2]];

    let mut m_inv = [[Complex::new(0.0, 0.0); 2]; 2];
    let det = m_in[0][0] * m_in[1][1] - m_in[0][1] * m_in[1][0];
    m_inv[0][0] = m_in[1][1] / det;
    m_inv[0][1] = -m_in[0][1] / det;
    m_inv[1][0] = -m_in[1][0] / det;
    m_inv[1][1] = m_in[0][0] / det;

    let mut m_temp = [[Complex::new(0.0, 0.0); 2]; 2];
    for j in 0..2 {
        for k in 0..2 {
            m_temp[j][k] = Complex::new(0.0, 0.0);
            for l in 0..2 {
                m_temp[j][k] += m[j][l] * m_out[l][k];
            }
        }
    }

    let mut m_total = [[Complex::new(0.0, 0.0); 2]; 2];
    for j in 0..2 {
        for k in 0..2 {
            m_total[j][k] = Complex::new(0.0, 0.0);
            for l in 0..2 {
                m_total[j][k] += m_inv[j][l] * m_temp[l][k];
            }
        }
    }

    let r_complex = m_total[1][0] / m_total[0][0];
    let t_complex = Complex::new(1.0, 0.0) / m_total[0][0];
    *r = r_complex.norm_sqr();
    *t = t_complex.norm_sqr() * (n_l * theta_l.cos()).re / (n_0 * theta.cos()).re;
}

fn main() {
    let start = Instant::now();
    let mut r = 0.0;
    let mut t = 0.0;

    // Define the layers manually
    let layers = vec![
        [Complex::new(1.5, 0.0), Complex::new(0.0, 0.0)], // Example layer 1
        [Complex::new(2.0, 0.0), Complex::new(100.0, 0.0)], // Example layer 2
        [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)], // Example layer 3 (substrate)
    ];

    solve_tmm(&mut r, &mut t, &layers, 500.0, 0.0);
    println!("Reflectance: {}", r);
    println!("Transmittance: {}", t);

    let elapsed = start.elapsed();
    println!("Execution Time: {:.6} seconds", elapsed.as_secs_f64());
}
