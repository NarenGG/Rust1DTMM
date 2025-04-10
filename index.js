import init, { initSync } from './pkg/Rust1DTMM.js';

async function run() {
    await init(); // Initialize the Wasm module
    console.log(initSync(5)); // Call the Rust function from JavaScript
}

run();
