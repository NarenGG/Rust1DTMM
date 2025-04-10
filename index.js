import init, { initSync } from './pkg/Rust1DTMM.js';

async function run() {
    await init(); // Initialize the Wasm module
    console.log(initSync()); // Call the Rust function from JavaScript
}

run();

export { init }; // Add a named export for `init`
export { solve_tmm } from './pkg/Rust1DTMM.js'; // Re-export `solve_tmm` from Rust1DTMM.js
export default init; // Keep the default export