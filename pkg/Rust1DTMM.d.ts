/* tslint:disable */
/* eslint-disable */
export function solve_tmm_js(layers: Float64Array, wavelength: number, theta: number): TMMResult;
export class TMMResult {
  private constructor();
  free(): void;
  reflectance: number;
  transmittance: number;
}

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
  readonly memory: WebAssembly.Memory;
  readonly __wbg_tmmresult_free: (a: number, b: number) => void;
  readonly __wbg_get_tmmresult_reflectance: (a: number) => number;
  readonly __wbg_set_tmmresult_reflectance: (a: number, b: number) => void;
  readonly __wbg_get_tmmresult_transmittance: (a: number) => number;
  readonly __wbg_set_tmmresult_transmittance: (a: number, b: number) => void;
  readonly solve_tmm_js: (a: number, b: number, c: number, d: number) => number;
  readonly __wbindgen_export_0: WebAssembly.Table;
  readonly __wbindgen_malloc: (a: number, b: number) => number;
  readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;
/**
* Instantiates the given `module`, which can either be bytes or
* a precompiled `WebAssembly.Module`.
*
* @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
*
* @returns {InitOutput}
*/
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
* If `module_or_path` is {RequestInfo} or {URL}, makes a request and
* for everything else, calls `WebAssembly.instantiate` directly.
*
* @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
*
* @returns {Promise<InitOutput>}
*/
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
