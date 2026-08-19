#include "book.h"
#include "lock.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define imin(a, b) (a < b ? a : b)

// Pascal launch configuration used by the GTX 1060 optimized version.
// The SM count is detected at runtime, so the same source also adapts to
// other Pascal GPUs when compiled for the appropriate architecture.
const int GTX1060_VECTOR_THREADS = 256;
const int GTX1060_VECTOR_BLOCKS_PER_SM = 8;

// Colored MATVEC/MTTVEC do not use atomics, but they remain the heaviest
// kernels.  Use the same 256-thread / 4-blocks-per-SM launch policy as the
// optimized atomic GTX 1060 implementation.
const int GTX1060_MATVEC_THREADS = 256;
const int GTX1060_MATVEC_BLOCKS_PER_SM = 4;

// Reduction block size must be a compile-time power of two because the
// kernels use a statically sized shared-memory reduction buffer.
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

/* Generator kolorów dla siatki HEX8:
   NAX, NAY, NAZ – liczba wezlów w kierunkach X, Y, Z
   ILE = (NAX-1)*(NAY-1)*(NAZ-1) – liczba elementów
   COLOR[e] – kolor elementu e (0..7)
*/

void generate_hex8_colors(int NAX, int NAY, int NAZ, int *COLOR)
{
//    int ILE = (NAX - 1) * (NAY - 1) * (NAZ - 1);
    int e = 0;
    int i, j, k;

    /* Zakladamy numeracje elementów:
       e = i + (NAX-1)*j + (NAX-1)*(NAY-1)*k
       i = 0..NAX-2, j = 0..NAY-2, k = 0..NAZ-2
    */

    for (k = 0; k < NAZ - 1; k++) {
        for (j = 0; j < NAY - 1; j++) {
            for (i = 0; i < NAX - 1; i++) {
                e = i + (NAX - 1) * j + (NAX - 1) * (NAY - 1) * k;

                /* 8-coloring:
                   kolor = (i mod 2) + 2*(j mod 2) + 4*(k mod 2)
                   daje kolory 0..7
                */
                COLOR[e] = (i & 1) + 2 * (j & 1) + 4 * (k & 1);
            }
        }
    }
}



/* Build compact element lists for the eight colors.
   COLOR_ELEMENTS[color_offsets[c] ... color_offsets[c+1]-1] contains only
   elements of color c.  This avoids scanning all ILE elements in every
   colored kernel launch. */
void build_color_element_lists(int ILE, const int *COLOR,
                               int *COLOR_ELEMENTS,
                               int color_offsets[9],
                               int color_counts[8])
{
    int cursor[8];

    for (int c = 0; c < 8; c++)
        color_counts[c] = 0;

    for (int e = 0; e < ILE; e++)
        color_counts[COLOR[e]]++;

    color_offsets[0] = 0;
    for (int c = 0; c < 8; c++)
        color_offsets[c + 1] = color_offsets[c] + color_counts[c];

    for (int c = 0; c < 8; c++)
        cursor[c] = color_offsets[c];

    for (int e = 0; e < ILE; e++) {
        int c = COLOR[e];
        COLOR_ELEMENTS[cursor[c]++] = e;
    }
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

__global__ void CZWARTA_cu(double *BK, double *BKNUM, double *BKDEN) {
  BK[0] = BKNUM[0] / BKDEN[0];
}

__global__ void CZWARTA2_cu(double *BKDEN, double *BKNUM) {
  BKDEN[0] = BKNUM[0];
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

__global__ void DRUGA_cu(Lock lock, int N, double *a, double *b, double *c) {
  __shared__ double cache[threadsPerBlock];
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int cacheIndex = threadIdx.x;

  double temp = 0;
  while (tid < N) {
    temp += a[tid] * b[tid];
    tid += blockDim.x * gridDim.x;
  }

  // Ustawienie wartoci w pamiàci podràcznej
  cache[cacheIndex] = temp;

  // Synchronizacja w¦tk¡w w tym bloku
  __syncthreads();

  // W przypadku redukcji threadsPerBlock musi byŠ potàg¦ 2,
  // ze wzglàdu na poni¬szy kod
  int i = blockDim.x / 2;
  while (i != 0) {
    if (cacheIndex < i)
      cache[cacheIndex] += cache[cacheIndex + i];
    __syncthreads();
    i /= 2;
  }

  if (cacheIndex == 0) {
    // Poczekanie na blokadà
    lock.lock();
    // Mamy blokadà, wiàc aktualizujemy i j¦ zwalniamy
    *c += cache[0];
    lock.unlock();
  }
}

__global__ void DNRM2_C_cu(Lock lock, int N, double *a, double *b) {
  __shared__ double cache[threadsPerBlock];
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int cacheIndex = threadIdx.x;

  double temp = 0;
  while (tid < N) {
    temp += pow(a[tid], 2);
    tid += blockDim.x * gridDim.x;
  }

  // Ustawienie wartoci w pamiàci podràcznej
  cache[cacheIndex] = temp;

  // Synchronizacja w¦tk¡w w tym bloku
  __syncthreads();

  // W przypadku redukcji threadsPerBlock musi byŠ potàg¦ 2,
  // ze wzglàdu na poni¬szy kod
  int i = blockDim.x / 2;
  while (i != 0) {
    if (cacheIndex < i)
      cache[cacheIndex] += cache[cacheIndex + i];
    __syncthreads();
    i /= 2;
  }

  if (cacheIndex == 0) {
    // Poczekanie na blokadà
    lock.lock();
    // Mamy blokadà, wiàc aktualizujemy i j¦ zwalniamy
    *b += cache[0];
    lock.unlock();
  }
}

__global__ void msolve(int N, double *R, double *Z, double *DIAG) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    Z[tid] = R[tid] / DIAG[tid];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void MATVEC_C_cu(Lock lock, double *XIN, double *XOUT, int *ISSp,
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

__global__ void MTTVEC_C_cu(Lock lock, double *XIN, double *XOUT, int *ISSp,
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

__global__ void MATVEC_color(
    const double *__restrict__ XIN,
    double *__restrict__ XOUT,
    const int *__restrict__ ISS,
    const double *__restrict__ BB,
    const int *__restrict__ COLOR_ELEMENTS,
    int color_offset, int color_count, int ILE)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    while (idx < color_count)
    {
        const int e = COLOR_ELEMENTS[color_offset + idx];
        int node[8];
        double xin[8];

        // Load connectivity and input values once per element.
        #pragma unroll
        for (int J = 0; J < 8; J++) {
            node[J] = ISS[J * ILE + e] - 1;
            xin[J] = XIN[node[J]];
        }

        // Coloring guarantees that elements in this launch do not share nodes,
        // so one accumulated store per output node is safe without atomicAdd.
        #pragma unroll
        for (int I = 0; I < 8; I++) {
            double acc = XOUT[node[I]];
            #pragma unroll
            for (int J = 0; J < 8; J++)
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
    int color_offset, int color_count, int ILE)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    while (idx < color_count)
    {
        const int e = COLOR_ELEMENTS[color_offset + idx];
        int node[8];
        double xin[8];

        #pragma unroll
        for (int J = 0; J < 8; J++) {
            node[J] = ISS[J * ILE + e] - 1;
            xin[J] = XIN[node[J]];
        }

        #pragma unroll
        for (int I = 0; I < 8; I++) {
            double acc = XOUT[node[I]];
            #pragma unroll
            for (int J = 0; J < 8; J++)
                acc += BB[I + 8 * J + 64 * e] * xin[J];
            XOUT[node[I]] = acc;
        }

        idx += blockDim.x * gridDim.x;
    }
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


extern "C" void kernel_dbcg_simple_(int *Np, double *B, double *X, double *TOLp,
                                    int *ITMAXp, int *ITERp, double *ERRp,
                                    double *BB, int *ISS, double *DIAG,
                                    int *NAXp, int *NAYp, int *NAZp,
									int *bc_dirichlet) {

  int N = *Np;
  double TOL = *TOLp;
  int ITMAX = *ITMAXp;
  int ITER = *ITERp;
  double ERR = *ERRp;
  int NAX = *NAXp;
  int NAY = *NAYp;
  int NAZ = *NAZp;

  int ILE = (NAX - 1) * (NAY - 1) * (NAZ - 1);
  int ILW = NAX * NAY * NAZ;

  int K;

  int *COLOR = (int *)malloc(ILE * sizeof(int));
  int *COLOR_ELEMENTS = (int *)malloc(ILE * sizeof(int));
  int color_offsets[9];
  int color_counts[8];

  double DNRM2, BNRM;
  //	  double AK, AKDEN, BK, BKDEN, BKNUM, BNRM, I;

  // CUDA part
  //****************************************************************************************
  // Detect the active GPU size at runtime, as in the GTX 1060 optimized version.
  int cuda_device = 0;
  cudaDeviceProp cuda_prop;
  cudaGetDevice(&cuda_device);
  cudaGetDeviceProperties(&cuda_prop, cuda_device);
  int sm_count = cuda_prop.multiProcessorCount;
  if (sm_count < 1) sm_count = 1;

  // Simple vector kernels.
  int Vthreads = GTX1060_VECTOR_THREADS;
  int Vblocks = imin(sm_count * GTX1060_VECTOR_BLOCKS_PER_SM,
                     (N + Vthreads - 1) / Vthreads);
  if (Vblocks < 1) Vblocks = 1;

  // Colored MATVEC/MTTVEC.  Each color gets its own grid sized to the compact
  // color list, capped at four blocks per SM.
  int Mthreads = GTX1060_MATVEC_THREADS;
  int color_blocks[8];

  // Dirichlet kernels.
  int Dthreads = GTX1060_VECTOR_THREADS;
  int Dblocks = imin(sm_count * GTX1060_VECTOR_BLOCKS_PER_SM,
                     (ILW + Dthreads - 1) / Dthreads);
  if (Dblocks < 1) Dblocks = 1;

  // Reduction kernels.
  int Rthreads = threadsPerBlock;
  int Rblocks = imin(sm_count * GTX1060_REDUCE_BLOCKS_PER_SM,
                     (N + Rthreads - 1) / Rthreads);
  if (Rblocks < 1) Rblocks = 1;


  double *P_cu, *PP_cu, *R_cu, *RR_cu, *Z_cu,
      *ZZ_cu; // declare GPU vector copies

  double *X_cu, *B_cu, *BB_cu, *DIAG_cu;

  int *ISS_cu;
  
  int *bc_dirichlet_cu;
 
  double *BNRM_cu, *DNRM2_cu;

  double *AK_cu, *AKDEN_cu, *BK_cu, *BKDEN_cu, *BKNUM_cu;

  Lock lock;

  cudaStream_t stream1, stream2;

  //*****************************************************************************************

  double *P, *PP, *R, *RR, *Z, *ZZ;

  P = (double *)calloc(N, sizeof(double));
  PP = (double *)calloc(N, sizeof(double));
  R = (double *)calloc(N, sizeof(double));
  RR = (double *)calloc(N, sizeof(double));
  Z = (double *)calloc(N, sizeof(double));
  ZZ = (double *)calloc(N, sizeof(double));
  /*
          P = (double *) malloc( N * sizeof(double) );
          PP= (double *) malloc( N * sizeof(double) );
          R = (double *) malloc( N * sizeof(double) );
          RR= (double *) malloc( N * sizeof(double) );
          Z = (double *) malloc( N * sizeof(double) );
          ZZ= (double *) malloc( N * sizeof(double) );
  */

  //     Cuda variables allocation on GPU
  //******************************************
  //     Allocate memory on GPU

  cudaMalloc((void **)&P_cu, sizeof(double) * N);
  cudaMalloc((void **)&R_cu, sizeof(double) * N);
  cudaMalloc((void **)&Z_cu, sizeof(double) * N);
  cudaMalloc((void **)&PP_cu, sizeof(double) * N);
  cudaMalloc((void **)&RR_cu, sizeof(double) * N);
  cudaMalloc((void **)&ZZ_cu, sizeof(double) * N);

  cudaMalloc((void **)&X_cu, sizeof(double) * ILW);
  cudaMalloc((void **)&B_cu, sizeof(double) * ILW);

  cudaMalloc((void **)&BB_cu, sizeof(double) * ILE * 64);
  cudaMalloc((void **)&ISS_cu, sizeof(int) * ILE * 8);
  cudaMalloc((void **)&DIAG_cu, sizeof(double) * ILW);
  
  cudaMalloc((void **)&bc_dirichlet_cu, sizeof(int) * ILW);
  

  cudaMalloc((void **)&BNRM_cu, sizeof(double));
  cudaMalloc((void **)&DNRM2_cu, sizeof(double));

  cudaMalloc((void **)&AK_cu, sizeof(double));
  cudaMalloc((void **)&AKDEN_cu, sizeof(double));
  cudaMalloc((void **)&BK_cu, sizeof(double));
  cudaMalloc((void **)&BKDEN_cu, sizeof(double));
  cudaMalloc((void **)&BKNUM_cu, sizeof(double));

  // copy vectors from CPU to GPU
  /*
        cudaMemcpy( P_cu,  P, sizeof(double) * N, cudaMemcpyHostToDevice );
        cudaMemcpy( R_cu,  R, sizeof(double) * N, cudaMemcpyHostToDevice );
        cudaMemcpy( Z_cu,  Z, sizeof(double) * N, cudaMemcpyHostToDevice );
        cudaMemcpy( PP_cu, PP, sizeof(double) * N, cudaMemcpyHostToDevice );
        cudaMemcpy( RR_cu, RR, sizeof(double) * N, cudaMemcpyHostToDevice );
        cudaMemcpy( ZZ_cu, ZZ, sizeof(double) * N, cudaMemcpyHostToDevice );
  */

  cudaMemcpy(X_cu, X, sizeof(double) * ILW, cudaMemcpyHostToDevice);
  cudaMemcpy(B_cu, B, sizeof(double) * ILW, cudaMemcpyHostToDevice);
  cudaMemcpy(BB_cu, BB, sizeof(double) * ILE * 64, cudaMemcpyHostToDevice);
  cudaMemcpy(ISS_cu, ISS, sizeof(int) * ILE * 8, cudaMemcpyHostToDevice);
  cudaMemcpy(DIAG_cu, DIAG, sizeof(double) * ILW, cudaMemcpyHostToDevice);
  
  cudaMemcpy(bc_dirichlet_cu, bc_dirichlet, sizeof(int) * ILW, cudaMemcpyHostToDevice);

  // seting memory value on GPU
  cudaMemset(P_cu, 0, sizeof(double) * N);
  cudaMemset(R_cu, 0, sizeof(double) * N);
  cudaMemset(Z_cu, 0, sizeof(double) * N);
  cudaMemset(PP_cu, 0, sizeof(double) * N);
  cudaMemset(RR_cu, 0, sizeof(double) * N);
  cudaMemset(ZZ_cu, 0, sizeof(double) * N);

  //******************************************

  cudaMemset(R_cu, 0, sizeof(double) * N);

  generate_hex8_colors(NAX, NAY, NAZ, COLOR);
  build_color_element_lists(ILE, COLOR, COLOR_ELEMENTS,
                            color_offsets, color_counts);

  // Precompute one launch size per compact color list.
  for (int col = 0; col < 8; col++) {
    color_blocks[col] = imin(sm_count * GTX1060_MATVEC_BLOCKS_PER_SM,
                             (color_counts[col] + Mthreads - 1) / Mthreads);
    if (color_blocks[col] < 1) color_blocks[col] = 1;
  }

  int *COLOR_ELEMENTS_cu;
  cudaMalloc(&COLOR_ELEMENTS_cu, ILE * sizeof(int));
  cudaMemcpy(COLOR_ELEMENTS_cu, COLOR_ELEMENTS, ILE * sizeof(int),
             cudaMemcpyHostToDevice);

  // Host coloring arrays are no longer needed after the compact list is copied.
  free(COLOR);
  free(COLOR_ELEMENTS);
  COLOR = NULL;
  COLOR_ELEMENTS = NULL;


  apply_dirichlet_border<<<Dblocks, Dthreads>>>(bc_dirichlet_cu, X_cu, ILW);
  
//  MATVEC_C_cu<<<blocks, Nth>>>(lock, X_cu, R_cu, ISS_cu, BB_cu, NAX, NAY, NAZ);

for (int col = 0; col < 8; col++)
{
    if (color_counts[col] == 0) continue;
    MATVEC_color<<<color_blocks[col], Mthreads>>>(
        X_cu,
        R_cu,
        ISS_cu,
        BB_cu,
        COLOR_ELEMENTS_cu,
        color_offsets[col],
        color_counts[col],
        ILE);
}




  
  apply_dirichlet<<<Dblocks, Dthreads>>>(bc_dirichlet_cu, R_cu, ILW);

  PIERWSZA_cu<<<Vblocks, Vthreads>>>(N, R_cu, B_cu, RR_cu);
  msolve<<<Vblocks, Vthreads>>>(N, R_cu, Z_cu, DIAG_cu);
  msolve<<<Vblocks, Vthreads>>>(N, RR_cu, ZZ_cu, DIAG_cu);

  cudaMemset(BNRM_cu, 0, sizeof(double));
  cudaMemset(DNRM2_cu, 0, sizeof(double));

  cudaStreamCreate(&stream1);
  cudaStreamCreate(&stream2);
  DNRM2_C_cu<<<Rblocks, Rthreads, 0, stream1>>>(lock, N, B_cu,
                                                             BNRM_cu);
  DNRM2_C_cu<<<Rblocks, Rthreads, 0, stream2>>>(lock, N, R_cu,
                                                             DNRM2_cu);
  cudaStreamSynchronize(stream1);
  cudaStreamSynchronize(stream2);
  cudaStreamDestroy(stream1);
  cudaStreamDestroy(stream2);

  cudaMemcpy(&DNRM2, DNRM2_cu, sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(&BNRM, BNRM_cu, sizeof(double), cudaMemcpyDeviceToHost);

  ERR = sqrt(DNRM2 / BNRM);

  if (ERR <= TOL)
    goto wyjscie;

  for (K = 1; K <= ITMAX; K++) {
    ITER = K;
    *ITERp = ITER;

    cudaMemset(BKNUM_cu, 0, sizeof(double));
    DRUGA_cu<<<Rblocks, Rthreads>>>(lock, N, Z_cu, RR_cu,
                                                 BKNUM_cu);

    if (ITER == 1) {
      TRZECIA_cu<<<Vblocks, Vthreads>>>(N, P_cu, Z_cu, PP_cu, ZZ_cu);
    } else {
      CZWARTA_cu<<<1, 1>>>(BK_cu, BKNUM_cu, BKDEN_cu);
      PIATA_cu<<<Vblocks, Vthreads>>>(N, BK_cu, P_cu, Z_cu, PP_cu, ZZ_cu);
    }

    CZWARTA2_cu<<<1, 1>>>(BKDEN_cu, BKNUM_cu);

    cudaMemset(Z_cu, 0, sizeof(double) * N);
	
	apply_dirichlet_border<<<Dblocks, Dthreads>>>(bc_dirichlet_cu, P_cu, ILW);
	
//    MATVEC_C_cu<<<blocks, Nth>>>(lock, P_cu, Z_cu, ISS_cu, BB_cu, NAX, NAY,
//                                 NAZ);

for (int col = 0; col < 8; col++)
{
    if (color_counts[col] == 0) continue;
    MATVEC_color<<<color_blocks[col], Mthreads>>>(
        P_cu,
        Z_cu,
        ISS_cu,
        BB_cu,
        COLOR_ELEMENTS_cu,
        color_offsets[col],
        color_counts[col],
        ILE);
}
								 
	apply_dirichlet<<<Dblocks, Dthreads>>>(bc_dirichlet_cu, Z_cu, ILW);							 

    cudaMemset(AKDEN_cu, 0, sizeof(double));
    DRUGA_cu<<<Rblocks, Rthreads>>>(lock, N, Z_cu, PP_cu,
                                                 AKDEN_cu);
    CZWARTA_cu<<<1, 1>>>(AK_cu, BKNUM_cu, AKDEN_cu);
    SZOSTA_cu<<<Vblocks, Vthreads>>>(N, AK_cu, X_cu, P_cu, R_cu, Z_cu);

    cudaMemset(ZZ_cu, 0, sizeof(double) * N);
	
	apply_dirichlet_border<<<Dblocks, Dthreads>>>(bc_dirichlet_cu, PP_cu, ILW);
	
//    MTTVEC_C_cu<<<blocks, Nth>>>(lock, PP_cu, ZZ_cu, ISS_cu, BB_cu, NAX, NAY,
//                                 NAZ);

for (int col = 0; col < 8; col++)
{
    if (color_counts[col] == 0) continue;
    MTTVEC_color<<<color_blocks[col], Mthreads>>>(
        PP_cu,
        ZZ_cu,
        ISS_cu,
        BB_cu,
        COLOR_ELEMENTS_cu,
        color_offsets[col],
        color_counts[col],
        ILE);
}

    apply_dirichlet<<<Dblocks, Dthreads>>>(bc_dirichlet_cu, ZZ_cu, ILW);

    SIODMA_cu<<<Vblocks, Vthreads>>>(N, AK_cu, RR_cu, ZZ_cu);

    cudaStreamCreate(&stream1);
    cudaStreamCreate(&stream2);

    msolve<<<Vblocks, Vthreads, 0, stream1>>>(N, R_cu, Z_cu, DIAG_cu);
    msolve<<<Vblocks, Vthreads, 0, stream2>>>(N, RR_cu, ZZ_cu, DIAG_cu);

    cudaStreamSynchronize(stream1);
    cudaStreamSynchronize(stream2);
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);

    cudaMemset(DNRM2_cu, 0, sizeof(double));
    DNRM2_C_cu<<<Rblocks, Rthreads>>>(lock, N, R_cu, DNRM2_cu);

    cudaMemcpy(&DNRM2, DNRM2_cu, sizeof(double), cudaMemcpyDeviceToHost);
    ERR = sqrt(DNRM2 / BNRM);
    *ERRp = ERR;

  //         printf(" ITER= %d ERR = %16.12e \n",ITER,ERR);
    if (ERR <= TOL)
      goto wyjscie;
  }

  //  printf("ERROR - stopping criterion not satisfied\n");

wyjscie:

  // copy solution to main function
  cudaMemcpy(X, X_cu, sizeof(double) * N, cudaMemcpyDeviceToHost);

  // free CPU memory

  free(P);
  free(PP);
  free(R);
  free(RR);
  free(Z);
  free(ZZ);

  // free GPU memory
  cudaFree(P_cu);
  cudaFree(R_cu);
  cudaFree(Z_cu);
  cudaFree(PP_cu);
  cudaFree(RR_cu);
  cudaFree(ZZ_cu);

  cudaFree(X_cu);
  cudaFree(B_cu);
  cudaFree(BB_cu);
  cudaFree(ISS_cu);
  cudaFree(DIAG_cu);
  
  cudaFree(bc_dirichlet_cu); 

  cudaFree(BNRM_cu);
  cudaFree(DNRM2_cu);

  cudaFree(AK_cu);
  cudaFree(AKDEN_cu);
  cudaFree(BK_cu);
  cudaFree(BKDEN_cu);
  cudaFree(BKNUM_cu);
  
  cudaFree(COLOR_ELEMENTS_cu);
  
  return;
}