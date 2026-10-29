# CUDA solver optimization experiments — Vanadis v2026.3.1

This directory contains **experimental CUDA implementations** of the
double-precision (FP64) DBCG solver used by Vanadis. The files are intended for
development, numerical comparison, and performance benchmarking; they are **not
part of the packaged release** and do not replace the default solver.

## Implementations

| Source file | Description | Validation status |
| --- | --- | --- |
| `cudatest_atomic8.cu` | Element-by-Element (EbE) HEX8 matrix-vector products with one thread per element and **8 FP64 atomic additions per element**, instead of 64. DBCG coefficients are calculated on the host. | Tested on a 2,800-node / 2,166-element case; results and iteration counts matched the earlier atomic implementation. |
| `cudatest_atomic8_gpu_scalars.cu` | Same atomic8 EbE kernels, with DBCG **alpha and beta calculated on the GPU**. The host still checks device-side breakdown/status after each coefficient and checks the residual norm for convergence. | GPU-scalar variant tested on the same case, with matching iteration counts and mass balance. The subsequent change to native FP64 `atomicAdd()` was reviewed statically but not rerun as a separate regression test. |

Both variants retain the original DBCG algorithm, Jacobi preconditioner,
convergence criterion, and Fortran entry point `kernel_dbcg_simple_`.

## Build and run

These are alternative implementations of the **same** CUDA entry point. Build
**only one** of them at a time, using the existing Vanadis Fortran sources and
build configuration. For an existing development Makefile that compiles
`cudatest.cu`, copy the desired variant to that filename in a separate build
directory (do not overwrite your canonical source without keeping a copy).

For the GTX 1060, compile with CUDA architecture `sm_61`, for example:

```sh
nvcc -O3 -arch=sm_61 -c cudatest.cu -o cudatest.o
```

Link the CUDA object with the existing Vanadis Fortran object(s), libraries,
and OpenMP options as specified by the development Makefile. A CUDA toolchain
and an NVIDIA GPU are required. These examples have been exercised with the
GTX 1060; compatibility and performance on other architectures should be
verified independently.

**Compatibility note:** native `atomicAdd(double*, double)` requires compute
capability 6.0 or newer. A retained CAS-based `atomicAddD` helper does **not**
automatically enable older GPUs if the kernels still call native `atomicAdd`.
Supporting pre-Pascal GPUs would require switching the relevant call sites
and building with an appropriate toolkit that supports those targets.

## Reference measurements

One test problem: 2,800 nodes, 2,166 HEX8 elements, 71 solver invocations,
1,593 DBCG iterations in total. Reported sums of solver times:

| Variant | Sum of solver times | Sum excluding first invocation |
| --- | ---: | ---: |
| Earlier atomic implementation (64 atomic updates/element) | 0.4744 s | 0.3777 s |
| `atomic8` | 0.4644 s | 0.3129 s |
| `atomic8_gpu_scalars` (before the native-atomic edit) | 0.4458 s | 0.3171 s |

The numerical results agreed: all 71 calls succeeded, with the same iteration
count for each call and effectively identical final mass balance. Timing
differences are **not conclusive performance claims**: GPU initialization,
first-call overhead, system load, and run-to-run variability matter. For
performance conclusions, repeat the measurements and separate startup time
from steady-state iterations.

## Benchmarking checklist

When adding a new variant, record:

1. GPU model, compute capability, CUDA toolkit, compiler flags, and CPU.
2. Mesh dimensions, number of HEX8 elements, and number of nodes.
3. Convergence tolerance, iteration counts, and final residual norms.
4. Total solver time, first-call time, and steady-state solver time.
5. Mass balance and differences relative to a reference solver.

The reductions use FP64 atomics, so floating-point accumulation order may
vary slightly between executions. Compare results using documented numerical
tolerances rather than requiring bitwise-identical output.

## Possible next experiments

- Move breakdown/status checks to the device and reduce host/device
  synchronization to the residual check, while preserving DBCG failure logic.
- Evaluate node-centric **gather** EbE matrix-vector multiplication without
  atomic writes, using a precomputed node-to-element adjacency.
- Investigate coalesced storage layouts for HEX8 element matrices and node
  indices.
- Fuse lightweight vector kernels or precompute the Jacobi inverse where
  numerical changes are acceptable and validated.

Keep these experiments in this directory. Promote an implementation into the
standard solver only after suitable numerical regression and performance tests.
