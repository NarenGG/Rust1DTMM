import init, { initSync } from './pkg/Rust1DTMM.js';

async function run() {
    await init(); // Initialize the Wasm module
    console.log(initSync()); // Call the Rust function from JavaScript
}

run();

export default init; // Add this line to provide a default export