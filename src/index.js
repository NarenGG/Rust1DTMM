import init, { greet } from './pkg/your_project_name.js';

async function run() {
    await init(); // Initialize the WebAssembly module
    console.log(greet("World")); // Call the exposed Rust function
}

run();
