import init, { initSync, solve_tmm_js } from './pkg/Rust1DTMM.js';

async function run() {
    await init(); // Initialize the Wasm module
    console.log(initSync()); // Call the Rust function from JavaScript
}

run();

export { init }; // Add a named export for `init`
export { solve_tmm_js as solve_tmm }; // Correctly re-export `solve_tmm_js` as `solve_tmm`
export default init; // Keep the default export