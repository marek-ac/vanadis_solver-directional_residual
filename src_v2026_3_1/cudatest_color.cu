#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define imin(a, b) ((a) < (b) ? (a) : (b))

// GTX 1060 / Pascal launch configuration.  The actual SM count is detected
// at runtime, so the same source adapts to both common GTX 1060 variants.
const int GTX1060_VECTOR_THREADS = 256;
const int GTX1060_VECTOR_BLOCKS_PER_SM = 8;

// Colored MATVEC/MTTVEC are the main optimization in this version: elements
// of one color do not share nodes, so EBE accumulation needs no atomics.
const int GTX1060_MATVEC_THREADS = 256;
const int GTX1060_MATVEC_BLOCKS_PER_SM = 4;

// Reduction block size must be a compile-time power of two.
const int threadsPerBlock = 512;
const int GTX1060_REDUCE_BLOCKS_PER_SM = 4;

__device__ double atomicAddD(double *address, double val) {
  unsigned long long int *address_as_ull =
      (unsigned long long int *)address;
  unsigned long long int old = *address_as_ull, assumed;

  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val +
                    __longlong_as_double(assumed)));
  } while (assumed != old);

  return __longlong_as_double(old);
}

// 8-coloring for a structured HEX8 mesh.
void generate_hex8_colors(int NAX, int NAY, int NAZ, int *COLOR) {
  int e, i, j, k;

  for (k = 0; k < NAZ - 1; ++k) {
    for (j = 0; j < NAY - 1; ++j) {
      for (i = 0; i < NAX - 1; ++i) {
        e = i + (NAX - 1) * j +
            (NAX - 1) * (NAY - 1) * k;
        COLOR[e] = (i & 1) + 2 * (j & 1) + 4 * (k & 1);
      }
    }
  }
}

// Compact element lists avoid scanning all elements for every color launch.
void build_color_element_lists(int ILE, const int *COLOR,
                               int *COLOR_ELEMENTS,
                               int color_offsets[9],
                               int color_counts[8]) {
  int cursor[8];

  for (int c = 0; c < 8; ++c)
    color_counts[c] = 0;

  for (int e = 0; e < ILE; ++e)
    color_counts[COLOR[e]]++;

  color_offsets[0] = 0;
  for (int c = 0; c < 8; ++c)
    color_offsets[c + 1] = color_offsets[c] + color_counts[c];

  for (int c = 0; c < 8; ++c)
    cursor[c] = color_offsets[c];

  for (int e = 0; e < ILE; ++e) {
    const int c = COLOR[e];
    COLOR_ELEMENTS[cursor[c]++] = e;
  }
}

__global__ void PIERWSZA_cu(int N, double *R, const double *B, double *RR) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    R[tid] = B[tid] - R[tid];
    RR[tid] = R[tid];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void PIATA_cu(int N, const double *BK, double *P,
                         const double *Z, double *PP, const double *ZZ) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    P[tid] = Z[tid] + BK[0] * P[tid];
    PP[tid] = ZZ[tid] + BK[0] * PP[tid];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void SIODMA_cu(int N, const double *AK, double *RR,
                          const double *ZZ) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    RR[tid] = RR[tid] - AK[0] * ZZ[tid];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void SZOSTA_cu(int N, const double *AK, double *X,
                          const double *P, double *R, const double *Z) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    X[tid] = X[tid] + AK[0] * P[tid];
    R[tid] = R[tid] - AK[0] * Z[tid];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void TRZECIA_cu(int N, double *P, const double *Z,
                           double *PP, const double *ZZ) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    P[tid] = Z[tid];
    PP[tid] = ZZ[tid];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void DRUGA_cu(int N, const double *a, const double *b,
                         double *c) {
  __shared__ double cache[threadsPerBlock];
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  const int cacheIndex = threadIdx.x;
  double temp = 0.0;

  while (tid < N) {
    temp += a[tid] * b[tid];
    tid += blockDim.x * gridDim.x;
  }

  cache[cacheIndex] = temp;
  __syncthreads();

  for (int i = blockDim.x / 2; i != 0; i /= 2) {
    if (cacheIndex < i)
      cache[cacheIndex] += cache[cacheIndex + i];
    __syncthreads();
  }

  if (cacheIndex == 0)
    atomicAddD(c, cache[0]);
}

__global__ void DNRM2_C_cu(int N, const double *a, double *b) {
  __shared__ double cache[threadsPerBlock];
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  const int cacheIndex = threadIdx.x;
  double temp = 0.0;

  while (tid < N) {
    temp += a[tid] * a[tid];
    tid += blockDim.x * gridDim.x;
  }

  cache[cacheIndex] = temp;
  __syncthreads();

  for (int i = blockDim.x / 2; i != 0; i /= 2) {
    if (cacheIndex < i)
      cache[cacheIndex] += cache[cacheIndex + i];
    __syncthreads();
  }

  if (cacheIndex == 0)
    atomicAddD(b, cache[0]);
}

__global__ void msolve(int N, const double *R, double *Z,
                       const double *DIAG) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    Z[tid] = R[tid] / DIAG[tid];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void MATVEC_color(
    const double *__restrict__ XIN,
    double *__restrict__ XOUT,
    const int *__restrict__ ISS,
    const double *__restrict__ BB,
    const int *__restrict__ COLOR_ELEMENTS,
    int color_offset, int color_count, int ILE) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;

  while (idx < color_count) {
    const int e = COLOR_ELEMENTS[color_offset + idx];
    int node[8];
    double xin[8];

#pragma unroll
    for (int J = 0; J < 8; ++J) {
      node[J] = ISS[J * ILE + e] - 1;
      xin[J] = XIN[node[J]];
    }

    // Same-color HEX8 elements do not share nodes.  The eight accumulated
    // stores are therefore race-free and need no atomicAdd.
#pragma unroll
    for (int I = 0; I < 8; ++I) {
      double acc = XOUT[node[I]];
#pragma unroll
      for (int J = 0; J < 8; ++J)
        acc += BB[J + 8 * I + 64 * e] * xin[J];
      XOUT[node[I]] = acc;
    }

    idx += blockDim.x * gridDim.x;
  }
}

__global__ void MTTVEC_color(
    const double *__restrict__ XIN,
    double *__restrict__ XOUT,
    const int *__restrict__ ISS,
    const double *__restrict__ BB,
    const int *__restrict__ COLOR_ELEMENTS,
    int color_offset, int color_count, int ILE) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;

  while (idx < color_count) {
    const int e = COLOR_ELEMENTS[color_offset + idx];
    int node[8];
    double xin[8];

#pragma unroll
    for (int J = 0; J < 8; ++J) {
      node[J] = ISS[J * ILE + e] - 1;
      xin[J] = XIN[node[J]];
    }

#pragma unroll
    for (int I = 0; I < 8; ++I) {
      double acc = XOUT[node[I]];
#pragma unroll
      for (int J = 0; J < 8; ++J)
        acc += BB[I + 8 * J + 64 * e] * xin[J];
      XOUT[node[I]] = acc;
    }

    idx += blockDim.x * gridDim.x;
  }
}

__global__ void apply_dirichlet_border(
    const int *__restrict__ bc_dirichlet,
    double *__restrict__ xin, int n) {
  for (int j = blockIdx.x * blockDim.x + threadIdx.x;
       j < n; j += blockDim.x * gridDim.x) {
    if (bc_dirichlet[j] == 1)
      xin[j] = 0.0;
  }
}

__global__ void apply_dirichlet(
    const int *__restrict__ bc_dirichlet,
    double *__restrict__ xout, int n) {
  for (int j = blockIdx.x * blockDim.x + threadIdx.x;
       j < n; j += blockDim.x * gridDim.x) {
    if (bc_dirichlet[j] == 1)
      xout[j] = 0.0;
  }
}

extern "C" void kernel_dbcg_simple_(int *Np, double *B, double *X,
                                    double *TOLp, int *ITMAXp,
                                    int *ITERp, double *ERRp,
                                    double *BB, int *ISS, double *DIAG,
                                    int *NAXp, int *NAYp, int *NAZp,
                                    int *bc_dirichlet, int *ISTATp) {
  // All declarations are intentionally before the first possible goto.
  // This keeps the error-cleanup path valid C++/CUDA code.
  const int N = *Np;
  const double TOL = *TOLp;
  const int ITMAX = *ITMAXp;
  const int NAX = *NAXp;
  const int NAY = *NAYp;
  const int NAZ = *NAZp;
  const int ILE = (NAX - 1) * (NAY - 1) * (NAZ - 1);
  const int ILW = NAX * NAY * NAZ;
  const double breakdown_tol = 1.0e-30;
  const double rhs_zero_tol2 = 1.0e-60;

  int status = 0;
  int ITER = 0;
  double ERR = 0.0;
  double DNRM2 = 0.0;
  double BNRM = 0.0;
  double BKNUM = 0.0;
  double BKDEN = 0.0;
  double BK = 0.0;
  double AKDEN = 0.0;
  double AK = 0.0;

  int cuda_device = 0;
  cudaDeviceProp cuda_prop;
  int sm_count = 1;
  int Vthreads = GTX1060_VECTOR_THREADS;
  int Vblocks = 1;
  int Mthreads = GTX1060_MATVEC_THREADS;
  int Dthreads = GTX1060_VECTOR_THREADS;
  int Dblocks = 1;
  int Rthreads = threadsPerBlock;
  int Rblocks = 1;
  int color_blocks[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int color_offsets[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  int color_counts[8] = {0, 0, 0, 0, 0, 0, 0, 0};

  int *COLOR = NULL;
  int *COLOR_ELEMENTS = NULL;
  int *COLOR_ELEMENTS_cu = NULL;
  int *ISS_cu = NULL;
  int *bc_dirichlet_cu = NULL;

  double *P_cu = NULL;
  double *PP_cu = NULL;
  double *R_cu = NULL;
  double *RR_cu = NULL;
  double *Z_cu = NULL;
  double *ZZ_cu = NULL;
  double *X_cu = NULL;
  double *B_cu = NULL;
  double *BB_cu = NULL;
  double *DIAG_cu = NULL;
  double *BNRM_cu = NULL;
  double *DNRM2_cu = NULL;
  double *AK_cu = NULL;
  double *BK_cu = NULL;
  double *AKDEN_cu = NULL;
  double *BKNUM_cu = NULL;

  *ITERp = 0;
  *ERRp = 0.0;
  *ISTATp = 0;

#define CUDA_CHECK(call) do {                                                   \
    cudaError_t _e = (call);                                                    \
    if (_e != cudaSuccess) {                                                    \
      fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__,        \
              cudaGetErrorString(_e));                                          \
      status = 5;                                                               \
      ERR = HUGE_VAL;                                                           \
      goto cleanup;                                                             \
    }                                                                           \
  } while (0)

#define CUDA_LAUNCH_CHECK() CUDA_CHECK(cudaPeekAtLastError())

  if (N <= 0 || NAX < 2 || NAY < 2 || NAZ < 2 || ILE <= 0 || ILW <= 0) {
    fprintf(stderr, "DBCG error: invalid mesh/system dimensions.\n");
    status = 5;
    ERR = HUGE_VAL;
    goto cleanup;
  }

  if (!isfinite(TOL) || TOL < 0.0 || ITMAX < 0) {
    fprintf(stderr, "DBCG error: invalid tolerance/iteration limit.\n");
    status = 5;
    ERR = HUGE_VAL;
    goto cleanup;
  }

  // Validate the host-side linear system before any GPU work.  A zero or
  // non-finite Jacobi diagonal would otherwise create NaN/Inf in msolve.
  for (int i = 0; i < N; ++i) {
    if (!isfinite(B[i]) || !isfinite(X[i])) {
      fprintf(stderr, "DBCG error: non-finite B/X at index %d\n", i + 1);
      status = 4;
      ERR = HUGE_VAL;
      goto cleanup;
    }
    if (!isfinite(DIAG[i]) || fabs(DIAG[i]) <= breakdown_tol) {
      fprintf(stderr,
              "DBCG error: invalid Jacobi diagonal at index %d: %.17e\n",
              i + 1, DIAG[i]);
      status = 2;
      ERR = HUGE_VAL;
      goto cleanup;
    }
  }

  COLOR = (int *)malloc((size_t)ILE * sizeof(int));
  COLOR_ELEMENTS = (int *)malloc((size_t)ILE * sizeof(int));
  if (COLOR == NULL || COLOR_ELEMENTS == NULL) {
    fprintf(stderr, "DBCG error: host allocation failed for coloring.\n");
    status = 5;
    ERR = HUGE_VAL;
    goto cleanup;
  }

  generate_hex8_colors(NAX, NAY, NAZ, COLOR);
  build_color_element_lists(ILE, COLOR, COLOR_ELEMENTS,
                            color_offsets, color_counts);
  if (color_offsets[8] != ILE) {
    fprintf(stderr, "DBCG error: invalid HEX8 color partition.\n");
    status = 5;
    ERR = HUGE_VAL;
    goto cleanup;
  }

  CUDA_CHECK(cudaGetDevice(&cuda_device));
  CUDA_CHECK(cudaGetDeviceProperties(&cuda_prop, cuda_device));
  sm_count = cuda_prop.multiProcessorCount;
  if (sm_count < 1)
    sm_count = 1;

  Vblocks = imin(sm_count * GTX1060_VECTOR_BLOCKS_PER_SM,
                  (N + Vthreads - 1) / Vthreads);
  if (Vblocks < 1)
    Vblocks = 1;

  Dblocks = imin(sm_count * GTX1060_VECTOR_BLOCKS_PER_SM,
                  (ILW + Dthreads - 1) / Dthreads);
  if (Dblocks < 1)
    Dblocks = 1;

  Rblocks = imin(sm_count * GTX1060_REDUCE_BLOCKS_PER_SM,
                  (N + Rthreads - 1) / Rthreads);
  if (Rblocks < 1)
    Rblocks = 1;

  for (int col = 0; col < 8; ++col) {
    if (color_counts[col] > 0) {
      color_blocks[col] = imin(sm_count * GTX1060_MATVEC_BLOCKS_PER_SM,
                               (color_counts[col] + Mthreads - 1) /
                               Mthreads);
      if (color_blocks[col] < 1)
        color_blocks[col] = 1;
    }
  }

  CUDA_CHECK(cudaMalloc((void **)&P_cu, sizeof(double) * (size_t)N));
  CUDA_CHECK(cudaMalloc((void **)&R_cu, sizeof(double) * (size_t)N));
  CUDA_CHECK(cudaMalloc((void **)&Z_cu, sizeof(double) * (size_t)N));
  CUDA_CHECK(cudaMalloc((void **)&PP_cu, sizeof(double) * (size_t)N));
  CUDA_CHECK(cudaMalloc((void **)&RR_cu, sizeof(double) * (size_t)N));
  CUDA_CHECK(cudaMalloc((void **)&ZZ_cu, sizeof(double) * (size_t)N));
  CUDA_CHECK(cudaMalloc((void **)&X_cu, sizeof(double) * (size_t)ILW));
  CUDA_CHECK(cudaMalloc((void **)&B_cu, sizeof(double) * (size_t)ILW));
  CUDA_CHECK(cudaMalloc((void **)&BB_cu,
                        sizeof(double) * (size_t)ILE * 64u));
  CUDA_CHECK(cudaMalloc((void **)&ISS_cu,
                        sizeof(int) * (size_t)ILE * 8u));
  CUDA_CHECK(cudaMalloc((void **)&DIAG_cu,
                        sizeof(double) * (size_t)ILW));
  CUDA_CHECK(cudaMalloc((void **)&bc_dirichlet_cu,
                        sizeof(int) * (size_t)ILW));
  CUDA_CHECK(cudaMalloc((void **)&COLOR_ELEMENTS_cu,
                        sizeof(int) * (size_t)ILE));
  CUDA_CHECK(cudaMalloc((void **)&BNRM_cu, sizeof(double)));
  CUDA_CHECK(cudaMalloc((void **)&DNRM2_cu, sizeof(double)));
  CUDA_CHECK(cudaMalloc((void **)&AK_cu, sizeof(double)));
  CUDA_CHECK(cudaMalloc((void **)&BK_cu, sizeof(double)));
  CUDA_CHECK(cudaMalloc((void **)&AKDEN_cu, sizeof(double)));
  CUDA_CHECK(cudaMalloc((void **)&BKNUM_cu, sizeof(double)));

  CUDA_CHECK(cudaMemcpy(X_cu, X, sizeof(double) * (size_t)ILW,
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(B_cu, B, sizeof(double) * (size_t)ILW,
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(BB_cu, BB,
                        sizeof(double) * (size_t)ILE * 64u,
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(ISS_cu, ISS,
                        sizeof(int) * (size_t)ILE * 8u,
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(DIAG_cu, DIAG,
                        sizeof(double) * (size_t)ILW,
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(bc_dirichlet_cu, bc_dirichlet,
                        sizeof(int) * (size_t)ILW,
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(COLOR_ELEMENTS_cu, COLOR_ELEMENTS,
                        sizeof(int) * (size_t)ILE,
                        cudaMemcpyHostToDevice));

  // Host coloring arrays are no longer needed after the compact list is on GPU.
  free(COLOR);
  free(COLOR_ELEMENTS);
  COLOR = NULL;
  COLOR_ELEMENTS = NULL;

  CUDA_CHECK(cudaMemset(P_cu, 0, sizeof(double) * (size_t)N));
  CUDA_CHECK(cudaMemset(R_cu, 0, sizeof(double) * (size_t)N));
  CUDA_CHECK(cudaMemset(Z_cu, 0, sizeof(double) * (size_t)N));
  CUDA_CHECK(cudaMemset(PP_cu, 0, sizeof(double) * (size_t)N));
  CUDA_CHECK(cudaMemset(RR_cu, 0, sizeof(double) * (size_t)N));
  CUDA_CHECK(cudaMemset(ZZ_cu, 0, sizeof(double) * (size_t)N));

  // Initial residual r = b - A*x.
  apply_dirichlet_border<<<Dblocks, Dthreads>>>(bc_dirichlet_cu,
                                                 X_cu, ILW);
  CUDA_LAUNCH_CHECK();

  for (int col = 0; col < 8; ++col) {
    if (color_counts[col] == 0)
      continue;
    MATVEC_color<<<color_blocks[col], Mthreads>>>(
        X_cu, R_cu, ISS_cu, BB_cu, COLOR_ELEMENTS_cu,
        color_offsets[col], color_counts[col], ILE);
    CUDA_LAUNCH_CHECK();
  }

  apply_dirichlet<<<Dblocks, Dthreads>>>(bc_dirichlet_cu, R_cu, ILW);
  CUDA_LAUNCH_CHECK();
  PIERWSZA_cu<<<Vblocks, Vthreads>>>(N, R_cu, B_cu, RR_cu);
  CUDA_LAUNCH_CHECK();
  msolve<<<Vblocks, Vthreads>>>(N, R_cu, Z_cu, DIAG_cu);
  CUDA_LAUNCH_CHECK();
  msolve<<<Vblocks, Vthreads>>>(N, RR_cu, ZZ_cu, DIAG_cu);
  CUDA_LAUNCH_CHECK();

  CUDA_CHECK(cudaMemset(BNRM_cu, 0, sizeof(double)));
  CUDA_CHECK(cudaMemset(DNRM2_cu, 0, sizeof(double)));
  DNRM2_C_cu<<<Rblocks, Rthreads>>>(N, B_cu, BNRM_cu);
  CUDA_LAUNCH_CHECK();
  DNRM2_C_cu<<<Rblocks, Rthreads>>>(N, R_cu, DNRM2_cu);
  CUDA_LAUNCH_CHECK();
  CUDA_CHECK(cudaMemcpy(&DNRM2, DNRM2_cu, sizeof(double),
                        cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(&BNRM, BNRM_cu, sizeof(double),
                        cudaMemcpyDeviceToHost));

  if (!isfinite(DNRM2) || !isfinite(BNRM) ||
      DNRM2 < 0.0 || BNRM < 0.0) {
    status = 4;
    ERR = HUGE_VAL;
    goto cleanup;
  }

  if (BNRM > rhs_zero_tol2)
    ERR = sqrt(DNRM2 / BNRM);
  else
    ERR = sqrt(DNRM2);  // absolute residual for zero RHS

  if (!isfinite(ERR)) {
    status = 4;
    ERR = HUGE_VAL;
    goto cleanup;
  }

  *ERRp = ERR;
  if (ERR <= TOL)
    goto copy_solution;

  for (int K = 1; K <= ITMAX; ++K) {
    ITER = K;
    *ITERp = ITER;

    CUDA_CHECK(cudaMemset(BKNUM_cu, 0, sizeof(double)));
    DRUGA_cu<<<Rblocks, Rthreads>>>(N, Z_cu, RR_cu, BKNUM_cu);
    CUDA_LAUNCH_CHECK();
    CUDA_CHECK(cudaMemcpy(&BKNUM, BKNUM_cu, sizeof(double),
                          cudaMemcpyDeviceToHost));

    if (!isfinite(BKNUM)) {
      status = 4;
      ERR = HUGE_VAL;
      goto cleanup;
    }

    if (K == 1) {
      TRZECIA_cu<<<Vblocks, Vthreads>>>(N, P_cu, Z_cu,
                                        PP_cu, ZZ_cu);
      CUDA_LAUNCH_CHECK();
    } else {
      if (fabs(BKDEN) <= breakdown_tol * fmax(1.0, fabs(BKNUM))) {
        status = 3;
        goto cleanup;
      }
      BK = BKNUM / BKDEN;
      if (!isfinite(BK)) {
        status = 4;
        ERR = HUGE_VAL;
        goto cleanup;
      }
      CUDA_CHECK(cudaMemcpy(BK_cu, &BK, sizeof(double),
                            cudaMemcpyHostToDevice));
      PIATA_cu<<<Vblocks, Vthreads>>>(N, BK_cu, P_cu, Z_cu,
                                     PP_cu, ZZ_cu);
      CUDA_LAUNCH_CHECK();
    }
    BKDEN = BKNUM;

    CUDA_CHECK(cudaMemset(Z_cu, 0, sizeof(double) * (size_t)N));
    apply_dirichlet_border<<<Dblocks, Dthreads>>>(bc_dirichlet_cu,
                                                   P_cu, ILW);
    CUDA_LAUNCH_CHECK();

    for (int col = 0; col < 8; ++col) {
      if (color_counts[col] == 0)
        continue;
      MATVEC_color<<<color_blocks[col], Mthreads>>>(
          P_cu, Z_cu, ISS_cu, BB_cu, COLOR_ELEMENTS_cu,
          color_offsets[col], color_counts[col], ILE);
      CUDA_LAUNCH_CHECK();
    }

    apply_dirichlet<<<Dblocks, Dthreads>>>(bc_dirichlet_cu, Z_cu, ILW);
    CUDA_LAUNCH_CHECK();

    CUDA_CHECK(cudaMemset(AKDEN_cu, 0, sizeof(double)));
    DRUGA_cu<<<Rblocks, Rthreads>>>(N, Z_cu, PP_cu, AKDEN_cu);
    CUDA_LAUNCH_CHECK();
    CUDA_CHECK(cudaMemcpy(&AKDEN, AKDEN_cu, sizeof(double),
                          cudaMemcpyDeviceToHost));

    if (!isfinite(AKDEN)) {
      status = 4;
      ERR = HUGE_VAL;
      goto cleanup;
    }
    if (fabs(AKDEN) <= breakdown_tol * fmax(1.0, fabs(BKNUM))) {
      status = 3;
      goto cleanup;
    }

    AK = BKNUM / AKDEN;
    if (!isfinite(AK)) {
      status = 4;
      ERR = HUGE_VAL;
      goto cleanup;
    }
    CUDA_CHECK(cudaMemcpy(AK_cu, &AK, sizeof(double),
                          cudaMemcpyHostToDevice));

    SZOSTA_cu<<<Vblocks, Vthreads>>>(N, AK_cu, X_cu,
                                     P_cu, R_cu, Z_cu);
    CUDA_LAUNCH_CHECK();

    CUDA_CHECK(cudaMemset(ZZ_cu, 0, sizeof(double) * (size_t)N));
    apply_dirichlet_border<<<Dblocks, Dthreads>>>(bc_dirichlet_cu,
                                                   PP_cu, ILW);
    CUDA_LAUNCH_CHECK();

    for (int col = 0; col < 8; ++col) {
      if (color_counts[col] == 0)
        continue;
      MTTVEC_color<<<color_blocks[col], Mthreads>>>(
          PP_cu, ZZ_cu, ISS_cu, BB_cu, COLOR_ELEMENTS_cu,
          color_offsets[col], color_counts[col], ILE);
      CUDA_LAUNCH_CHECK();
    }

    apply_dirichlet<<<Dblocks, Dthreads>>>(bc_dirichlet_cu, ZZ_cu, ILW);
    CUDA_LAUNCH_CHECK();
    SIODMA_cu<<<Vblocks, Vthreads>>>(N, AK_cu, RR_cu, ZZ_cu);
    CUDA_LAUNCH_CHECK();

    msolve<<<Vblocks, Vthreads>>>(N, R_cu, Z_cu, DIAG_cu);
    CUDA_LAUNCH_CHECK();
    msolve<<<Vblocks, Vthreads>>>(N, RR_cu, ZZ_cu, DIAG_cu);
    CUDA_LAUNCH_CHECK();

    CUDA_CHECK(cudaMemset(DNRM2_cu, 0, sizeof(double)));
    DNRM2_C_cu<<<Rblocks, Rthreads>>>(N, R_cu, DNRM2_cu);
    CUDA_LAUNCH_CHECK();
    CUDA_CHECK(cudaMemcpy(&DNRM2, DNRM2_cu, sizeof(double),
                          cudaMemcpyDeviceToHost));

    if (!isfinite(DNRM2) || DNRM2 < 0.0) {
      status = 4;
      ERR = HUGE_VAL;
      goto cleanup;
    }

    if (BNRM > rhs_zero_tol2)
      ERR = sqrt(DNRM2 / BNRM);
    else
      ERR = sqrt(DNRM2);

    if (!isfinite(ERR)) {
      status = 4;
      ERR = HUGE_VAL;
      goto cleanup;
    }

    *ERRp = ERR;
    if (ERR <= TOL)
      goto copy_solution;
  }

  status = 1;

copy_solution:
  CUDA_CHECK(cudaMemcpy(X, X_cu, sizeof(double) * (size_t)N,
                        cudaMemcpyDeviceToHost));

cleanup:
  *ITERp = ITER;
  *ERRp = ERR;
  *ISTATp = status;

  if (COLOR != NULL)
    free(COLOR);
  if (COLOR_ELEMENTS != NULL)
    free(COLOR_ELEMENTS);

  if (P_cu != NULL) cudaFree(P_cu);
  if (R_cu != NULL) cudaFree(R_cu);
  if (Z_cu != NULL) cudaFree(Z_cu);
  if (PP_cu != NULL) cudaFree(PP_cu);
  if (RR_cu != NULL) cudaFree(RR_cu);
  if (ZZ_cu != NULL) cudaFree(ZZ_cu);
  if (X_cu != NULL) cudaFree(X_cu);
  if (B_cu != NULL) cudaFree(B_cu);
  if (BB_cu != NULL) cudaFree(BB_cu);
  if (ISS_cu != NULL) cudaFree(ISS_cu);
  if (DIAG_cu != NULL) cudaFree(DIAG_cu);
  if (bc_dirichlet_cu != NULL) cudaFree(bc_dirichlet_cu);
  if (COLOR_ELEMENTS_cu != NULL) cudaFree(COLOR_ELEMENTS_cu);
  if (BNRM_cu != NULL) cudaFree(BNRM_cu);
  if (DNRM2_cu != NULL) cudaFree(DNRM2_cu);
  if (AK_cu != NULL) cudaFree(AK_cu);
  if (BK_cu != NULL) cudaFree(BK_cu);
  if (AKDEN_cu != NULL) cudaFree(AKDEN_cu);
  if (BKNUM_cu != NULL) cudaFree(BKNUM_cu);

#undef CUDA_CHECK
#undef CUDA_LAUNCH_CHECK
}
