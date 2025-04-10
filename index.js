import init, { solve_tmm_js } from './pkg/Rust1DTMM.js';

async function run() {
    await init(); // Initialize the Wasm module
    console.log(solve_tmm_js(5)); // Call the Rust function from JavaScript
}

run();
