#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "book.h"
#include "lock.h"

#define imin(a,b) (a<b?a:b)


const int threadsPerBlock = 128;
const int blocksPerGrid   = 16;

/*

$ gcc -c -O3 bdbcg.c
$ gfortran -c -O3 VANADIS.FOR
$ gfortran -o  vanadis.out VANADIS.o bdbcg.o
$ ./vanadis.out

*/


__device__ double atomicAddD(double* address, double val)
{
   unsigned long long int* address_as_ull =
                             (unsigned long long int*)address;
   unsigned long long int old = *address_as_ull, assumed;

   do {
       assumed = old;
       old = atomicCAS(address_as_ull, assumed,
                       __double_as_longlong(val +
                              __longlong_as_double(assumed)));

   // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
   } while (assumed != old);

   return __longlong_as_double(old);
}


__global__ void PIERWSZA_cu(int N, double *R, double *B, double *RR)
{
   int tid = threadIdx.x + blockIdx.x * blockDim.x;
   while (tid < N) {
     R[tid] = B[tid] - R[tid];
     RR[tid] = R[tid];
     tid += blockDim.x * gridDim.x;
   }
}

__global__ void PIATA_cu(int N, double *BK, double *P, double *Z, double *PP, double *ZZ)
{
   int tid = threadIdx.x + blockIdx.x * blockDim.x;
   while (tid < N) {
     P[tid]  = Z[tid]  + BK[0]*P[tid];
     PP[tid] = ZZ[tid] + BK[0]*PP[tid];
     tid += blockDim.x * gridDim.x;
   }
}


__global__ void SIODMA_cu(int N, double *AK, double *RR, double *ZZ)      
{
   int tid = threadIdx.x + blockIdx.x * blockDim.x;
   while (tid < N) {
     RR[tid]  = RR[tid] - AK[0]*ZZ[tid];
     tid += blockDim.x * gridDim.x;
   }
}

__global__ void SZOSTA_cu(int N, double *AK, double *X, double *P, double *R, double*Z)
{
   int tid = threadIdx.x + blockIdx.x * blockDim.x;
   while (tid < N) {
     X[tid] = X[tid] + AK[0]*P[tid];
     R[tid] = R[tid] - AK[0]*Z[tid];
     tid += blockDim.x * gridDim.x;
   }
}

__global__ void CZWARTA_cu(double *BK, double *BKNUM, double *BKDEN)
{
                BK[0]=BKNUM[0]/BKDEN[0];
}


__global__ void CZWARTA2_cu( double *BKDEN, double *BKNUM)
{
                BKDEN[0]=BKNUM[0];
}


__global__ void TRZECIA_cu(int N, double *P, double *Z, double *PP, double *ZZ)
{
   int tid = threadIdx.x + blockIdx.x * blockDim.x;
   while (tid < N) {
     P[tid]  = Z[tid];
     PP[tid] = ZZ[tid];
     tid += blockDim.x * gridDim.x;
   }
}


__global__ void DRUGA_cu( Lock lock, int N, double *a, double *b ,double *c) {
    __shared__ double cache[threadsPerBlock];
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int cacheIndex = threadIdx.x;

    double   temp = 0;
    while (tid < N) {
        temp += a[tid]*b[tid];
        tid += blockDim.x * gridDim.x;
    }
    
    // Ustawienie wartoci w pamiàci podràcznej
    cache[cacheIndex] = temp;
    
    // Synchronizacja w¦tk¡w w tym bloku
    __syncthreads();

    // W przypadku redukcji threadsPerBlock musi byŠ potàg¦ 2,
    // ze wzglàdu na poni¬szy kod
    int i = blockDim.x/2;
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




__global__ void DNRM2_C_cu( Lock lock, int N, double *a, double *b ) {
    __shared__ double cache[threadsPerBlock];
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int cacheIndex = threadIdx.x;

    double   temp = 0;
    while (tid < N) {
        temp += pow(a[tid],2);
        tid += blockDim.x * gridDim.x;
    }
    
    // Ustawienie wartoci w pamiàci podràcznej
    cache[cacheIndex] = temp;
    
    // Synchronizacja w¦tk¡w w tym bloku
    __syncthreads();

    // W przypadku redukcji threadsPerBlock musi byŠ potàg¦ 2,
    // ze wzglàdu na poni¬szy kod
    int i = blockDim.x/2;
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



// simple kernel function that adds two vectors
__global__ void vect_add(double *a, double *b, int N)
{
   int idx = threadIdx.x;
   if (idx<N) a[idx] = a[idx] + b[idx];
}



__global__ void msolve(int N, double *R, double *Z, double *DIAG)
{
   int tid = threadIdx.x + blockIdx.x * blockDim.x;
   while (tid < N) {
     Z[tid] = R[tid] / DIAG[tid];
     tid += blockDim.x * gridDim.x;
   }
   
}


// function called from main fortran program
extern "C" void kernel_msolve_(int *Np, double *R, double *Z, double *DIAG)
{
   double  *a_d, *b_d, *c_d;  // declare GPU vector copies
   
   int blocks = 128;     // uses 128 block of
   int Nth = 1024;       // N threads on GPU

   int N = * Np;
   
      

   // Allocate memory on GPU
   cudaMalloc( (void **)&a_d, sizeof(double) * N );
   cudaMalloc( (void **)&b_d, sizeof(double) * N );
   cudaMalloc( (void **)&c_d, sizeof(double) * N );

   // copy vectors from CPU to GPU
   cudaMemcpy( a_d, R, sizeof(double) * N, cudaMemcpyHostToDevice );
   cudaMemcpy( b_d, Z, sizeof(double) * N, cudaMemcpyHostToDevice );
   cudaMemcpy( c_d, DIAG, sizeof(double) * N, cudaMemcpyHostToDevice );


   // call function on GPU
   msolve<<< blocks, Nth >>>(N, a_d, b_d, c_d);

   // copy vectors back from GPU to CPU
   cudaMemcpy( R, a_d, sizeof(double) * N, cudaMemcpyDeviceToHost );
   cudaMemcpy( Z, b_d, sizeof(double) * N, cudaMemcpyDeviceToHost );
   cudaMemcpy( DIAG, c_d, sizeof(double) * N, cudaMemcpyDeviceToHost );


   // free GPU memory
   cudaFree(a_d);
   cudaFree(b_d);
   cudaFree(c_d);

   return;
}




      void MSOLVE_C(int N, double *R, double *Z, double* DIAG)
      {
	  int I;
      
		  for (I=0; I < N ; I++) 
                      Z[I] = R[I]/DIAG[I];
      return;		
      }
	  
      void MTSOLVE_C(int N, double *RR, double *ZZ, double* DIAG)
      {
	  int I;
      
		  for (I=0; I < N ; I++)
			ZZ[I] = RR[I]/DIAG[I];
		  
      return;		
      }
	  

      void DNRM2_C(int N, double *DX, double * DNRM)
      {
        int I;

        double DNRM2 = 0;

        for (I=0; I < N ; I++) 
         {
           DNRM2 = DNRM2 + pow(DX[I],2);
   //        DNRM2 = DNRM2 + DX[I] * DX[I];

       }
	DNRM2 = sqrt(DNRM2);
       *DNRM = DNRM2;

        return;
      }

	  
__global__ void MATVEC_C_cu(Lock lock, double *XIN, double* XOUT,int* ISSp, double* BB, int NAX, int NAY, int NAZ)
      {
      
	  int ILOSC_EL = (NAX-1) * (NAY-1) * (NAZ-1);
          int I, J ;
              

  	  int tid = threadIdx.x + blockIdx.x * blockDim.x;
      
      
              while (tid < ILOSC_EL)

		  {
		    for (I=0; I < 8; I++)
		       for (J=0; J < 8; J++)
                          {
                  //           Fp[0] = XOUT[ISSp[I*ILOSC_EL+tid]-1];
                             atomicAddD(& XOUT[ISSp[I*ILOSC_EL+tid]-1]  ,  BB[J+8*I+64*tid]*XIN[ISSp[J*ILOSC_EL+tid]-1]);
                          }
                    tid += blockDim.x * gridDim.x;

		  }
	  
	  return;	  

      }
	  

__global__ void MATVEC_C_cuf(Lock lock, double *XIN, double* XOUT,int* ISSp, double* BB, int NAX, int NAY, int NAZ)
      {
      
	  int ILOSC_EL = (NAX-1) * (NAY-1) * (NAZ-1);
          int I, J ;

         __shared__  double temp_v[256][8];
         memset( temp_v,  0, sizeof(double)*256*8);

         __shared__  int temp_p[256][8];

         __syncthreads();
              
  	  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  
      
              while (tid < ILOSC_EL)

		  {
		       for (I=0; I < 8; I++)                       
                       {
                          temp_p[threadIdx.x][I] = ISSp[I*ILOSC_EL+tid]-1;

      		          for (J=0; J < 8; J++)
                          {
                             
                             atomicAddD(&temp_v[threadIdx.x][I], BB[J+8*I+64*tid] * XIN[ISSp[J*ILOSC_EL+tid]-1]);
//                             atomicAddD(& XOUT[ISSp[I*ILOSC_EL+tid]-1]  ,  BB[J+8*I+64*tid]*XIN[ISSp[J*ILOSC_EL+tid]-1]);

                          }
                         }
                    tid += blockDim.x * gridDim.x;

		  }

        __syncthreads();
	  

          for (I=0; I < 8; I++)
             atomicAddD(&XOUT[temp_p[threadIdx.x][I]], temp_v[threadIdx.x][I]);

	  
          return;	  

      }
	


  
__global__ void MTTVEC_C_cu(Lock lock, double *XIN, double* XOUT,int* ISSp, double* BB, int NAX, int NAY, int NAZ)
      {
      
	  int ILOSC_EL = (NAX-1) * (NAY-1) * (NAZ-1);
          int I, J ;
              
  	  int tid = threadIdx.x + blockIdx.x * blockDim.x;
      
      
              while (tid < ILOSC_EL)

		  {
		    for (I=0; I < 8; I++)
		       for (J=0; J < 8; J++)
                          atomicAddD(&XOUT[ISSp[I*ILOSC_EL+tid]-1], BB[I+8*J+64*tid]*XIN[ISSp[J*ILOSC_EL+tid]-1]);

                    tid += blockDim.x * gridDim.x;

		  }
	  
	  return;	  
      }


void MATVEC_C(double *XIN, double* XOUT,int* ISSp, double* BB, int NAX, int NAY, int NAZ)
      {
      

	  int I_OT_ELEM[8];
	  int ILOSC_EL = (NAX-1) * (NAY-1) * (NAZ-1);
	  int ILOSC_WEZLOW = NAX * NAY * NAZ;
          int I, J, IK, JE;
	  
          // int ISS[8][ILOSC_EL];
          // memcpy (ISS, ISSp, sizeof(ISS));

		  memset(XOUT, 0, ILOSC_WEZLOW *sizeof(double));
		  IK = 0;
		  
		  for (JE=1; JE <= ILOSC_EL; JE++)
		  {
		  
			for (J=0; J < 8; J++)
			//  I_OT_ELEM[J]=ISS[J][JE-1];
                           I_OT_ELEM[J]=ISSp[J*ILOSC_EL+JE-1];

			for (I=0; I < 8; I++)
			{
			  for (J=0; J < 8; J++)
			  {
				XOUT[I_OT_ELEM[I]-1]=XOUT[I_OT_ELEM[I]-1] +BB[IK]*XIN[I_OT_ELEM[J]-1];                                
				IK++;
			  }
			}
               

		  }
	  
	  return;
	  }
	  
	  
	  void MTTVEC_C(double *XIN, double* XOUT,int* ISSp, double* BB, int NAX, int NAY, int NAZ)
      {
      
	  int I_OT_ELEM[8];
	  int ILOSC_EL = (NAX-1) * (NAY-1) * (NAZ-1);
	  int ILOSC_WEZLOW = NAX * NAY * NAZ;
          int I, J, IK3, JE;
	  double B3[8][8];

//          int ISS[8][ILOSC_EL];
//          memcpy (ISS, ISSp, sizeof(ISS));

	  
		  memset(XOUT, 0, ILOSC_WEZLOW *sizeof(double));
		  
		  for (JE=0; JE < ILOSC_EL; JE++)
		  {
			for (J=0; J < 8; J++)
                      //     I_OT_ELEM[J]=ISS[J][JE];
			  I_OT_ELEM[J]=ISSp[J*ILOSC_EL+JE];

			IK3=0;
			for (I=0; I < 8; I++)
			{
                          IK3++;
 			 for (J=0; J < 8; J++)
		//	    B3[I][J]=BB[J+8*(IK3-1)+64*(JE-1)];
                            B3[I][J]=BB[J+8*(IK3-1)+64*(JE)];			  
			}
			
			for (I=0; I < 8; I++)
			{
			  for (J=0; J < 8; J++)
			  {
				XOUT[I_OT_ELEM[I]-1]=XOUT[I_OT_ELEM[I]-1] +B3[J][I]*XIN[I_OT_ELEM[J]-1];
			  }
			}
		  }
	  
	  return;
	  }



extern "C"  void kernel_dbcg_simple_(int *Np, double *B, double *X, double *TOLp, int *ITMAXp, int *ITERp,
	                      double *ERRp, double *BB, int *ISS, double *DIAG, int *NAXp, int *NAYp,int *NAZp)
{
	 


	  int N = *Np;
	  double TOL = *TOLp;
	  int ITMAX = *ITMAXp;
	  int ITER = *ITERp;
	  double ERR = *ERRp;
	  int NAX = *NAXp;
	  int NAY = *NAYp;
	  int NAZ = *NAZp;

          int ILE =  (NAX-1) * (NAY-1) * (NAZ-1); 
          int ILW =  NAX * NAY * NAZ;
        
        
	  int K;	  
	  double DNRM2, BNRM;
//	  double AK, AKDEN, BK, BKDEN, BKNUM, BNRM, I;

// CUDA part        
//****************************************************************************************
//      const int threadsPerBlock = 256;
//      const int blocksPerGrid   = imin( 32, (N+threadsPerBlock-1) / threadsPerBlock );

        int blocks = 256;     // uses 128 block of
        int Nth = 1024;       // N threads on GPU

        double  *P_cu, *PP_cu, *R_cu, *RR_cu, *Z_cu, *ZZ_cu;  // declare GPU vector copies

        double  *X_cu, *B_cu, *BB_cu, *DIAG_cu;
      
        int *ISS_cu;
      
        double *BNRM_cu, *DNRM2_cu;
      
        double *AK_cu, *AKDEN_cu, *BK_cu, *BKDEN_cu, *BKNUM_cu;
            
        Lock    lock;

        cudaStream_t stream1, stream2; 
 
        int GPU_CPU = 0;

     
//*****************************************************************************************

//        double P[N], PP[N], R[N], RR[N], Z[N], ZZ[N];


        double *P, *PP, *R, *RR, *Z, *ZZ;    
   

   
	P = (double *) calloc( N , sizeof(double) );
	PP= (double *) calloc( N , sizeof(double) );
	R = (double *) calloc( N , sizeof(double) );
	RR= (double *) calloc( N , sizeof(double) );
	Z = (double *) calloc( N , sizeof(double) );
	ZZ= (double *) calloc( N , sizeof(double) );
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

      cudaMalloc( (void **)&P_cu, sizeof(double) * N );
      cudaMalloc( (void **)&R_cu, sizeof(double) * N );
      cudaMalloc( (void **)&Z_cu, sizeof(double) * N );
      cudaMalloc( (void **)&PP_cu, sizeof(double) * N );
      cudaMalloc( (void **)&RR_cu, sizeof(double) * N );
      cudaMalloc( (void **)&ZZ_cu, sizeof(double) * N );

      cudaMalloc( (void **)&X_cu, sizeof(double) * ILW );
      cudaMalloc( (void **)&B_cu, sizeof(double) * ILW );

      cudaMalloc( (void **)&BB_cu, sizeof(double) * ILE * 64 );
      cudaMalloc( (void **)&ISS_cu, sizeof(int) * ILE * 8 );
      cudaMalloc( (void **)&DIAG_cu, sizeof(double) * ILW );
      
      cudaMalloc( (void **)&BNRM_cu, sizeof(double) );
      cudaMalloc( (void **)&DNRM2_cu, sizeof(double) );
            
      cudaMalloc( (void **)&AK_cu, sizeof(double) );
      cudaMalloc( (void **)&AKDEN_cu, sizeof(double) );
      cudaMalloc( (void **)&BK_cu, sizeof(double) );
      cudaMalloc( (void **)&BKDEN_cu, sizeof(double) );
      cudaMalloc( (void **)&BKNUM_cu, sizeof(double) );      



      
// copy vectors from CPU to GPU
/*
      cudaMemcpy( P_cu,  P, sizeof(double) * N, cudaMemcpyHostToDevice );
      cudaMemcpy( R_cu,  R, sizeof(double) * N, cudaMemcpyHostToDevice );
      cudaMemcpy( Z_cu,  Z, sizeof(double) * N, cudaMemcpyHostToDevice );
      cudaMemcpy( PP_cu, PP, sizeof(double) * N, cudaMemcpyHostToDevice );
      cudaMemcpy( RR_cu, RR, sizeof(double) * N, cudaMemcpyHostToDevice );
      cudaMemcpy( ZZ_cu, ZZ, sizeof(double) * N, cudaMemcpyHostToDevice );
*/



      cudaMemcpy( X_cu, X, sizeof(double) * ILW, cudaMemcpyHostToDevice );
      cudaMemcpy( B_cu, B, sizeof(double) * ILW, cudaMemcpyHostToDevice );
      cudaMemcpy( BB_cu, BB, sizeof(double) * ILE * 64, cudaMemcpyHostToDevice );
      cudaMemcpy( ISS_cu, ISS, sizeof(int) * ILE * 8, cudaMemcpyHostToDevice );
      cudaMemcpy( DIAG_cu, DIAG, sizeof(double) * ILW, cudaMemcpyHostToDevice );


// seting memory value on GPU
      cudaMemset( P_cu,   0, sizeof(double) * N);
      cudaMemset( R_cu,   0, sizeof(double) * N);
      cudaMemset( Z_cu,   0, sizeof(double) * N);
      cudaMemset( PP_cu,  0, sizeof(double) * N);
      cudaMemset( RR_cu,  0, sizeof(double) * N);
      cudaMemset( ZZ_cu,  0, sizeof(double) * N);
      
//******************************************        
        


      if (GPU_CPU != 1)
      {
        cudaMemset( R_cu,   0, sizeof(double) * N);   
        MATVEC_C_cu<<<blocks, Nth>>>(lock, X_cu, R_cu, ISS_cu, BB_cu, NAX, NAY, NAZ);
      }
      else
      {
         cudaMemcpy( X, X_cu, sizeof(double) * N, cudaMemcpyDeviceToHost );
         MATVEC_C(X, R, ISS, BB, NAX, NAY, NAZ);
         cudaMemcpy( R_cu, R, sizeof(double) * N, cudaMemcpyHostToDevice );
      }

      PIERWSZA_cu<<<blocks, Nth>>>(N, R_cu, B_cu, RR_cu);
      msolve<<< blocks, Nth >>>(N,R_cu,Z_cu,DIAG_cu);
      msolve<<< blocks, Nth >>>(N,RR_cu,ZZ_cu,DIAG_cu);




// test ******
// copy solution to main function
//      int I;
//      cudaMemcpy( ZZ, ZZ_cu, sizeof(double) * N, cudaMemcpyDeviceToHost );
//      for (I=0; I < N; I++){
//          if (ZZ[I] != 0) printf("***i = %d  ZZ[i] = %15.6e \n",I,ZZ[I]);
//      }
// test******

      cudaMemset( BNRM_cu,  0, sizeof(double));
      cudaMemset( DNRM2_cu, 0, sizeof(double));   

      cudaStreamCreate(&stream1);
      cudaStreamCreate(&stream2);
      DNRM2_C_cu<<<blocksPerGrid,threadsPerBlock,0,stream1>>>( lock, N, B_cu, BNRM_cu);
      DNRM2_C_cu<<<blocksPerGrid,threadsPerBlock,0,stream2>>>( lock, N, R_cu, DNRM2_cu);
      cudaStreamSynchronize(stream1);
      cudaStreamSynchronize(stream2);  	
      cudaStreamDestroy(stream1);
      cudaStreamDestroy(stream2);

      cudaMemcpy( &DNRM2, DNRM2_cu, sizeof(double) , cudaMemcpyDeviceToHost );
      cudaMemcpy( &BNRM,  BNRM_cu,  sizeof(double) , cudaMemcpyDeviceToHost );

      ERR = sqrt(DNRM2/BNRM);


      if (ERR <= TOL) goto wyjscie;
	  
      for (K = 1 ; K <= ITMAX; K++)
      {
         ITER  = K;
        *ITERp = ITER;

         cudaMemset( BKNUM_cu, 0, sizeof(double));
         DRUGA_cu<<<blocksPerGrid,threadsPerBlock>>>(lock, N, Z_cu, RR_cu, BKNUM_cu);

	 if(ITER == 1)
         {
            TRZECIA_cu<<<blocks, Nth>>>(N, P_cu, Z_cu, PP_cu, ZZ_cu);
	 }  
         else
	 {
            CZWARTA_cu<<<1,1>>>(BK_cu, BKNUM_cu, BKDEN_cu);
            PIATA_cu<<<blocks, Nth>>>(N, BK_cu, P_cu, Z_cu, PP_cu, ZZ_cu);
         }
         
         CZWARTA2_cu<<<1,1>>>(BKDEN_cu, BKNUM_cu);


         if (GPU_CPU != 1)
         {
           cudaMemset( Z_cu,   0, sizeof(double) * N); 
  	   MATVEC_C_cu<<<blocks, Nth>>>(lock,P_cu,Z_cu,ISS_cu,BB_cu,NAX,NAY,NAZ);
//           MATVEC_C_cuf<<<256, 256>>>(lock,P_cu,Z_cu,ISS_cu,BB_cu,NAX,NAY,NAZ);
         }
         else
         {
           cudaMemcpy( P, P_cu, sizeof(double) * N, cudaMemcpyDeviceToHost );
           MATVEC_C(P,Z,ISS,BB,NAX,NAY,NAZ);
           cudaMemcpy( Z_cu, Z, sizeof(double) * N, cudaMemcpyHostToDevice );
         }

         cudaMemset( AKDEN_cu, 0, sizeof(double));
         DRUGA_cu<<<blocksPerGrid,threadsPerBlock>>>(lock, N, Z_cu, PP_cu, AKDEN_cu);
         CZWARTA_cu<<<1,1>>>(AK_cu, BKNUM_cu, AKDEN_cu);
         SZOSTA_cu<<<blocks, Nth>>>(N, AK_cu, X_cu, P_cu, R_cu, Z_cu);         


         if (GPU_CPU != 1)
         {
	   cudaMemset( ZZ_cu,   0, sizeof(double) * N);	
           MTTVEC_C_cu<<<blocks, Nth>>>(lock,PP_cu,ZZ_cu,ISS_cu,BB_cu,NAX,NAY,NAZ);

         }
         else
         {
           cudaMemcpy( PP, PP_cu, sizeof(double) * N, cudaMemcpyDeviceToHost );
           MTTVEC_C(PP,ZZ,ISS,BB,NAX,NAY,NAZ);
           cudaMemcpy( ZZ_cu, ZZ, sizeof(double) * N, cudaMemcpyHostToDevice );
         }

         SIODMA_cu<<<blocks, Nth>>>(N, AK_cu, RR_cu, ZZ_cu);


         cudaStreamCreate(&stream1);
         cudaStreamCreate(&stream2);
      
         msolve<<< blocks, Nth, 0, stream1 >>>(N,R_cu,Z_cu,DIAG_cu);
         msolve<<< blocks, Nth, 0, stream2 >>>(N,RR_cu,ZZ_cu,DIAG_cu);

         cudaStreamSynchronize(stream1);
         cudaStreamSynchronize(stream2);  	
         cudaStreamDestroy(stream1);
         cudaStreamDestroy(stream2);

         cudaMemset( DNRM2_cu, 0, sizeof(double));
         DNRM2_C_cu<<<blocksPerGrid,threadsPerBlock>>>( lock, N, R_cu, DNRM2_cu);
         cudaMemcpy( &DNRM2,  DNRM2_cu,  sizeof(double) , cudaMemcpyDeviceToHost );
         ERR = sqrt(DNRM2/BNRM); 
        *ERRp = ERR;

//         printf(" ITER= %d ERR = %16.12e \n",ITER,ERR);
         if (ERR <= TOL) goto wyjscie;
       }
	  
	  
      printf("PROBLEM Z ROZWIAZANIEM UKLADU ROWNAN\n");

      wyjscie:
 
 // copy solution to main function
      cudaMemcpy( X, X_cu, sizeof(double) * N, cudaMemcpyDeviceToHost );

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
      
        cudaFree(BNRM_cu);
        cudaFree(DNRM2_cu);
      
        cudaFree(AK_cu);
        cudaFree(AKDEN_cu);
        cudaFree(BK_cu);
        cudaFree(BKDEN_cu);
        cudaFree(BKNUM_cu);           
      return;
      }




