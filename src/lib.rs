use wasm_bindgen::prelude::*;
use wasm_bindgen_rayon::init_thread_pool;
use num::complex::Complex;
use std::f64::consts::PI;
// ... your other imports ...

#[wasm_bindgen(start)]
pub async fn run() -> Result<(), JsValue> {
    // The number of threads to use. Typically, you'll want to match the number of cores.
    // num_cpus::get() gives the number of logical cores, but you might limit this value.
    init_thread_pool(num_cpus::get() as u32).await?;
    Ok(())
}

// Your existing `solve_tmm_js` function that uses Rayon can then use into_par_iter().
#[wasm_bindgen]
pub async fn solve_tmm_js(
    layers: Vec<f64>,
    wavelength: f64,
    theta: f64,
) -> Result<TMMResult, JsValue> {
    // Now Rayon’s parallel iterators will work since the thread pool is initialized.
    let num_layers = layers.len() / 3;
    let parsed_layers: Vec<[Complex<f64>; 2]> = (0..num_layers)
        .into_par_iter() // This is now allowed with Rayon in WASM.
        .map(|i| {
            let real_part = layers[i * 3];
            let imag_part = layers[i * 3 + 1];
            let thickness = layers[i * 3 + 2];
            [
                Complex::new(real_part, imag_part),
                Complex::new(thickness, 0.0),
            ]
        })
        .collect();

    let mut r = 0.0;
    let mut t = 0.0;

    solve_tmm(&mut r, &mut t, &parsed_layers, num_layers, wavelength, theta);

    Ok(TMMResult {
        reflectance: r,
        transmittance: t,
    })
}
// Transfer Matrix Method function (unchanged)
fn solve_tmm(
    r: &mut f64,
    t: &mut f64,
    layers: &Vec<[Complex<f64>; 2]>,
    num_layers: usize,
    wavelength: f64,
    theta: f64,
) {
    let k_0 = (2.0 * PI) / wavelength;

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

// Transfer Matrix helper function (unchanged)
fn transfer_matrix(
    result: &mut [[Complex<f64>; 2]; 2],
    k_0: f64,
    n: Complex<f64>,
    d: f64,
    theta: f64,
) {
    let k_z = k_0 * n * theta.cos();
    let q_1 = cos_complex(k_z * d);
    let q_2: Complex<f64> = Complex::<f64>::i() * sin_complex(k_z * d);
    let n_cos_th = n * theta.cos();

    result[0][0] = q_1;
    result[0][1] = q_2 / n_cos_th;
    result[1][0] = n_cos_th * q_2;
    result[1][1] = q_1;
}

// Complex sine and cosine functions (unchanged)
fn sin_complex(v: Complex<f64>) -> Complex<f64> {
    let a = v.re % (2.0 * PI);
    let b = v.im % (2.0 * PI);
    Complex::new(a.sin() * b.cosh(), a.cos() * b.sinh())
}

fn cos_complex(v: Complex<f64>) -> Complex<f64> {
    let a = v.re % (2.0 * PI);
    let b = v.im % (2.0 * PI);
    Complex::new(a.cos() * b.cosh(), -a.sin() * b.sinh())
}