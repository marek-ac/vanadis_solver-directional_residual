#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define imin(a, b) (a < b ? a : b)

// GeForce GTX 1060 (Pascal, compute capability 6.1) launch configuration.
// The actual SM count is detected at runtime, so this also adapts to both
// common GTX 1060 variants (6 GB and 3 GB) without hard-coding 10 or 9 SMs.
// All affected kernels use grid-stride loops.
const int GTX1060_VECTOR_THREADS = 256;
const int GTX1060_VECTOR_BLOCKS_PER_SM = 8;

// MATVEC/MTTVEC are heavy per thread: each element performs 64 FP64
// multiply/atomic-add contributions.  A 256-thread block and four blocks
// per SM is a conservative starting point for GP106/CC 6.1: enough warps
// to hide latency without creating an unnecessarily large atomic workload.
const int GTX1060_MATVEC_THREADS = 256;
const int GTX1060_MATVEC_BLOCKS_PER_SM = 4;

// Reduction kernels use a static shared-memory array, therefore the block
// size is a compile-time constant and must remain a power of two.
// Four 512-thread blocks per SM can expose up to 2048 threads/SM on CC 6.1,
// useful for the GTX 1060's relatively slow FP64 path, while still using
// fewer final lock acquisitions than the original fixed 64-block launch.
const int threadsPerBlock = 512;
const int GTX1060_REDUCE_BLOCKS_PER_SM = 4;

__device__ double atomicAddD(double *address, double val) {
  unsigned long long int *address_as_ull = (unsigned long long int *)address;
  unsigned long long int old = *address_as_ull, assumed;

  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val + __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN !=
    // NaN)
  } while (assumed != old);

  return __longlong_as_double(old);
}

__global__ void PIERWSZA_cu(int N, double *R, double *B, double *RR) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    R[tid] = B[tid] - R[tid];
    RR[tid] = R[tid];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void PIATA_cu(int N, double *BK, double *P, double *Z, double *PP,
                         double *ZZ) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    P[tid] = Z[tid] + BK[0] * P[tid];
    PP[tid] = ZZ[tid] + BK[0] * PP[tid];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void SIODMA_cu(int N, double *AK, double *RR, double *ZZ) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    RR[tid] = RR[tid] - AK[0] * ZZ[tid];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void SZOSTA_cu(int N, double *AK, double *X, double *P, double *R,
                          double *Z) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    X[tid] = X[tid] + AK[0] * P[tid];
    R[tid] = R[tid] - AK[0] * Z[tid];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void TRZECIA_cu(int N, double *P, double *Z, double *PP,
                           double *ZZ) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    P[tid] = Z[tid];
    PP[tid] = ZZ[tid];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void DRUGA_cu(int N, double *a, double *b, double *c) {
  __shared__ double cache[threadsPerBlock];
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int cacheIndex = threadIdx.x;

  double temp = 0;
  while (tid < N) {
    temp += a[tid] * b[tid];
    tid += blockDim.x * gridDim.x;
  }

  // Store the partial sum in shared memory
  cache[cacheIndex] = temp;

  // Synchronize threads in this block
  __syncthreads();

  // For this reduction threadsPerBlock must be a power of two,
  // as required by the reduction loop below
  int i = blockDim.x / 2;
  while (i != 0) {
    if (cacheIndex < i)
      cache[cacheIndex] += cache[cacheIndex + i];
    __syncthreads();
    i /= 2;
  }

  if (cacheIndex == 0) {
    atomicAddD(c, cache[0]);
  }
}

__global__ void DNRM2_C_cu(int N, double *a, double *b) {
  __shared__ double cache[threadsPerBlock];
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int cacheIndex = threadIdx.x;

  double temp = 0;
  while (tid < N) {
    temp += a[tid] * a[tid];
    tid += blockDim.x * gridDim.x;
  }

  // Store the partial sum in shared memory
  cache[cacheIndex] = temp;

  // Synchronize threads in this block
  __syncthreads();

  // For this reduction threadsPerBlock must be a power of two,
  // as required by the reduction loop below
  int i = blockDim.x / 2;
  while (i != 0) {
    if (cacheIndex < i)
      cache[cacheIndex] += cache[cacheIndex + i];
    __syncthreads();
    i /= 2;
  }

  if (cacheIndex == 0) {
    atomicAddD(b, cache[0]);
  }
}

__global__ void msolve(int N, double *R, double *Z, double *DIAG) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    Z[tid] = R[tid] / DIAG[tid];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void MATVEC_C_cu(double *XIN, double *XOUT, int *ISSp,
                            double *BB, int NAX, int NAY, int NAZ) {

  int ILOSC_EL = (NAX - 1) * (NAY - 1) * (NAZ - 1);
  int I, J;

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < ILOSC_EL)

  {
    for (I = 0; I < 8; I++)
      for (J = 0; J < 8; J++) {
        //           Fp[0] = XOUT[ISSp[I*ILOSC_EL+tid]-1];
        atomicAdd(&XOUT[ISSp[I * ILOSC_EL + tid] - 1],
                  BB[J + 8 * I + 64 * tid] * XIN[ISSp[J * ILOSC_EL + tid] - 1]);
      }
    tid += blockDim.x * gridDim.x;
  }

  return;
}

__global__ void MTTVEC_C_cu(double *XIN, double *XOUT, int *ISSp,
                            double *BB, int NAX, int NAY, int NAZ) {

  int ILOSC_EL = (NAX - 1) * (NAY - 1) * (NAZ - 1);
  int I, J;

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < ILOSC_EL)

  {
    for (I = 0; I < 8; I++)
      for (J = 0; J < 8; J++)
        atomicAdd(&XOUT[ISSp[I * ILOSC_EL + tid] - 1],
                  BB[I + 8 * J + 64 * tid] * XIN[ISSp[J * ILOSC_EL + tid] - 1]);

    tid += blockDim.x * gridDim.x;
  }

  return;
}

__global__ void apply_dirichlet_border(
    const int* __restrict__ bc_dirichlet,
    double* __restrict__ xin,
    int n)
{
    for (int j = blockIdx.x * blockDim.x + threadIdx.x;
         j < n;
         j += blockDim.x * gridDim.x)
    {
        if (bc_dirichlet[j] == 1)
            xin[j] = 0.0;
    }
}


__global__ void apply_dirichlet(
    const int* __restrict__ bc_dirichlet,
    double* __restrict__ xout,
    int n)
{
    for (int j = blockIdx.x * blockDim.x + threadIdx.x;
         j < n;
         j += blockDim.x * gridDim.x)
    {
        if (bc_dirichlet[j] == 1)
            xout[j] = 0.0;
    }
}


__global__ void apply_dirichlet_gstride(
    const int* __restrict__ bc_dirichlet,
    double* __restrict__ xout,
    int n)
{
    for (int j = blockIdx.x * blockDim.x + threadIdx.x;
         j < n;
         j += blockDim.x * gridDim.x)
    {
        if (bc_dirichlet[j] == 1)
            xout[j] = 0.0;
    }
}



__global__ void apply_dirichlet_branchless(
    const int* __restrict__ bc_dirichlet,
    double* __restrict__ xout,
    int n)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j < n) {
        int mask = (bc_dirichlet[j] == 1);
        xout[j] *= (1 - mask);   // jesli mask=1 ? xout=0
    }
}

__global__ void apply_dirichlet_ldg(
    const int* bc_dirichlet,
    double* xout,
    int n)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j < n) {
        int bc = __ldg(&bc_dirichlet[j]);
        if (bc == 1)
            xout[j] = 0.0;
    }
}


extern "C" void kernel_dbcg_simple_(int *Np, double *B, double *X,
                                    double *TOLp, int *ITMAXp,
                                    int *ITERp, double *ERRp,
                                    double *BB, int *ISS, double *DIAG,
                                    int *NAXp, int *NAYp, int *NAZp,
                                    int *bc_dirichlet, int *ISTATp) {

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
  double DNRM2 = 0.0, BNRM = 0.0;
  double BKNUM = 0.0, BKDEN = 0.0, BK = 0.0;
  double AKDEN = 0.0, AK = 0.0;

  double *P_cu = NULL, *PP_cu = NULL, *R_cu = NULL;
  double *RR_cu = NULL, *Z_cu = NULL, *ZZ_cu = NULL;
  double *X_cu = NULL, *B_cu = NULL, *BB_cu = NULL, *DIAG_cu = NULL;
  double *BNRM_cu = NULL, *DNRM2_cu = NULL;
  double *AK_cu = NULL, *BK_cu = NULL, *AKDEN_cu = NULL;
  double *BKNUM_cu = NULL;
  int *ISS_cu = NULL, *bc_dirichlet_cu = NULL;

  // All function-scope locals that may be crossed by a forward goto
  // must be declared before the first possible `goto cleanup` in C++.
  int cuda_device = 0;
  cudaDeviceProp cuda_prop;
  int sm_count = 1;
  int Vthreads = GTX1060_VECTOR_THREADS;
  int Vblocks = 1;
  int Mthreads = GTX1060_MATVEC_THREADS;
  int Mblocks = 1;
  int Dthreads = GTX1060_VECTOR_THREADS;
  int Dblocks = 1;
  int Rthreads = threadsPerBlock;
  int Rblocks = 1;

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

#define CUDA_LAUNCH_CHECK() do {                                                \
    cudaError_t _e = cudaGetLastError();                                        \
    if (_e != cudaSuccess) {                                                    \
      fprintf(stderr, "CUDA kernel launch error at %s:%d: %s\n",              \
              __FILE__, __LINE__, cudaGetErrorString(_e));                      \
      status = 5;                                                               \
      ERR = HUGE_VAL;                                                           \
      goto cleanup;                                                             \
    }                                                                           \
  } while (0)

  // Validate the host-side linear system before any GPU work.
  for (int i = 0; i < N; ++i) {
    if (!isfinite(B[i]) || !isfinite(X[i])) {
      fprintf(stderr, "DBCG error: non-finite B/X at index %d\n", i + 1);
      status = 4;
      ERR = HUGE_VAL;
      goto cleanup;
    }
    if (!isfinite(DIAG[i]) || fabs(DIAG[i]) <= breakdown_tol) {
      fprintf(stderr, "DBCG error: invalid Jacobi diagonal at index %d: %.17e\n",
              i + 1, DIAG[i]);
      status = 2;
      ERR = HUGE_VAL;
      goto cleanup;
    }
  }

  CUDA_CHECK(cudaGetDevice(&cuda_device));
  CUDA_CHECK(cudaGetDeviceProperties(&cuda_prop, cuda_device));
  sm_count = cuda_prop.multiProcessorCount;
  if (sm_count < 1) sm_count = 1;

  Vblocks = imin(sm_count * GTX1060_VECTOR_BLOCKS_PER_SM,
                 (N + Vthreads - 1) / Vthreads);
  if (Vblocks < 1) Vblocks = 1;

  Mblocks = imin(sm_count * GTX1060_MATVEC_BLOCKS_PER_SM,
                 (ILE + Mthreads - 1) / Mthreads);
  if (Mblocks < 1) Mblocks = 1;

  Dblocks = imin(sm_count * GTX1060_VECTOR_BLOCKS_PER_SM,
                 (ILW + Dthreads - 1) / Dthreads);
  if (Dblocks < 1) Dblocks = 1;

  Rblocks = imin(sm_count * GTX1060_REDUCE_BLOCKS_PER_SM,
                 (N + Rthreads - 1) / Rthreads);
  if (Rblocks < 1) Rblocks = 1;

  CUDA_CHECK(cudaMalloc((void **)&P_cu, sizeof(double) * N));
  CUDA_CHECK(cudaMalloc((void **)&R_cu, sizeof(double) * N));
  CUDA_CHECK(cudaMalloc((void **)&Z_cu, sizeof(double) * N));
  CUDA_CHECK(cudaMalloc((void **)&PP_cu, sizeof(double) * N));
  CUDA_CHECK(cudaMalloc((void **)&RR_cu, sizeof(double) * N));
  CUDA_CHECK(cudaMalloc((void **)&ZZ_cu, sizeof(double) * N));
  CUDA_CHECK(cudaMalloc((void **)&X_cu, sizeof(double) * ILW));
  CUDA_CHECK(cudaMalloc((void **)&B_cu, sizeof(double) * ILW));
  CUDA_CHECK(cudaMalloc((void **)&BB_cu, sizeof(double) * ILE * 64));
  CUDA_CHECK(cudaMalloc((void **)&ISS_cu, sizeof(int) * ILE * 8));
  CUDA_CHECK(cudaMalloc((void **)&DIAG_cu, sizeof(double) * ILW));
  CUDA_CHECK(cudaMalloc((void **)&bc_dirichlet_cu, sizeof(int) * ILW));
  CUDA_CHECK(cudaMalloc((void **)&BNRM_cu, sizeof(double)));
  CUDA_CHECK(cudaMalloc((void **)&DNRM2_cu, sizeof(double)));
  CUDA_CHECK(cudaMalloc((void **)&AK_cu, sizeof(double)));
  CUDA_CHECK(cudaMalloc((void **)&BK_cu, sizeof(double)));
  CUDA_CHECK(cudaMalloc((void **)&AKDEN_cu, sizeof(double)));
  CUDA_CHECK(cudaMalloc((void **)&BKNUM_cu, sizeof(double)));

  CUDA_CHECK(cudaMemcpy(X_cu, X, sizeof(double) * ILW,
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(B_cu, B, sizeof(double) * ILW,
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(BB_cu, BB, sizeof(double) * ILE * 64,
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(ISS_cu, ISS, sizeof(int) * ILE * 8,
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(DIAG_cu, DIAG, sizeof(double) * ILW,
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(bc_dirichlet_cu, bc_dirichlet,
                        sizeof(int) * ILW, cudaMemcpyHostToDevice));

  CUDA_CHECK(cudaMemset(P_cu, 0, sizeof(double) * N));
  CUDA_CHECK(cudaMemset(R_cu, 0, sizeof(double) * N));
  CUDA_CHECK(cudaMemset(Z_cu, 0, sizeof(double) * N));
  CUDA_CHECK(cudaMemset(PP_cu, 0, sizeof(double) * N));
  CUDA_CHECK(cudaMemset(RR_cu, 0, sizeof(double) * N));
  CUDA_CHECK(cudaMemset(ZZ_cu, 0, sizeof(double) * N));

  apply_dirichlet_border<<<Dblocks, Dthreads>>>(bc_dirichlet_cu, X_cu,
                                                 ILW);
  CUDA_LAUNCH_CHECK();
  MATVEC_C_cu<<<Mblocks, Mthreads>>>(X_cu, R_cu, ISS_cu, BB_cu,
                                     NAX, NAY, NAZ);
  CUDA_LAUNCH_CHECK();
  apply_dirichlet<<<Dblocks, Dthreads>>>(bc_dirichlet_cu, R_cu, ILW);
  CUDA_LAUNCH_CHECK();
  PIERWSZA_cu<<<Vblocks, Vthreads>>>(N, R_cu, B_cu, RR_cu);
  CUDA_LAUNCH_CHECK();
  msolve<<<Vblocks, Vthreads>>>(N, R_cu, Z_cu, DIAG_cu);
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

  if (!isfinite(DNRM2) || !isfinite(BNRM) || DNRM2 < 0.0 || BNRM < 0.0) {
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
      TRZECIA_cu<<<Vblocks, Vthreads>>>(N, P_cu, Z_cu, PP_cu, ZZ_cu);
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

    CUDA_CHECK(cudaMemset(Z_cu, 0, sizeof(double) * N));
    apply_dirichlet_border<<<Dblocks, Dthreads>>>(bc_dirichlet_cu, P_cu,
                                                   ILW);
    CUDA_LAUNCH_CHECK();
    MATVEC_C_cu<<<Mblocks, Mthreads>>>(P_cu, Z_cu, ISS_cu, BB_cu,
                                       NAX, NAY, NAZ);
    CUDA_LAUNCH_CHECK();
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

    SZOSTA_cu<<<Vblocks, Vthreads>>>(N, AK_cu, X_cu, P_cu, R_cu, Z_cu);
    CUDA_LAUNCH_CHECK();

    CUDA_CHECK(cudaMemset(ZZ_cu, 0, sizeof(double) * N));
    apply_dirichlet_border<<<Dblocks, Dthreads>>>(bc_dirichlet_cu, PP_cu,
                                                   ILW);
    CUDA_LAUNCH_CHECK();
    MTTVEC_C_cu<<<Mblocks, Mthreads>>>(PP_cu, ZZ_cu, ISS_cu, BB_cu,
                                       NAX, NAY, NAZ);
    CUDA_LAUNCH_CHECK();
    apply_dirichlet<<<Dblocks, Dthreads>>>(bc_dirichlet_cu, ZZ_cu, ILW);
    CUDA_LAUNCH_CHECK();
    SIODMA_cu<<<Vblocks, Vthreads>>>(N, AK_cu, RR_cu, ZZ_cu);
    CUDA_LAUNCH_CHECK();

    msolve<<<Vblocks, Vthreads>>>(N, R_cu, Z_cu, DIAG_cu);
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
  goto copy_solution;

copy_solution:
  CUDA_CHECK(cudaMemcpy(X, X_cu, sizeof(double) * N,
                        cudaMemcpyDeviceToHost));

cleanup:
  *ITERp = ITER;
  *ERRp = ERR;
  *ISTATp = status;

  if (P_cu) cudaFree(P_cu);
  if (R_cu) cudaFree(R_cu);
  if (Z_cu) cudaFree(Z_cu);
  if (PP_cu) cudaFree(PP_cu);
  if (RR_cu) cudaFree(RR_cu);
  if (ZZ_cu) cudaFree(ZZ_cu);
  if (X_cu) cudaFree(X_cu);
  if (B_cu) cudaFree(B_cu);
  if (BB_cu) cudaFree(BB_cu);
  if (ISS_cu) cudaFree(ISS_cu);
  if (DIAG_cu) cudaFree(DIAG_cu);
  if (bc_dirichlet_cu) cudaFree(bc_dirichlet_cu);
  if (BNRM_cu) cudaFree(BNRM_cu);
  if (DNRM2_cu) cudaFree(DNRM2_cu);
  if (AK_cu) cudaFree(AK_cu);
  if (BK_cu) cudaFree(BK_cu);
  if (AKDEN_cu) cudaFree(AKDEN_cu);
  if (BKNUM_cu) cudaFree(BKNUM_cu);

#undef CUDA_CHECK
#undef CUDA_LAUNCH_CHECK
  return;
}

