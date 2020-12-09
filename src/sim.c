/*
 * sim.c - similarity (half-)matrix
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>         // timing
#include <R_ext/Applic.h> // blas routines

#include "etc.h"

void
xTx(
  int M,
  int N,
  double *X,
  int Cstride,
  double *C
  )
{
    // use optimized blas routine to do X^T X
    char *transN = "N", *transT= "T";
    double one = 1.0;
    double zero = 0.0;
    double *Ctmp = (double*)nalloc((size_t)N*N,sizeof(double));

    F77_CALL(dgemm)(transT, transN, &N, &N, &M, &one,
     X, &M, X, &M, &zero, Ctmp, &N); 

    for(size_t i = 0; i < N; i++ )
      {
      double *Ci = C + (i*(i+1)/2) * Cstride, *Ctmpi = Ctmp + i*N;
      for(int j = 0; j <= i; j++ )
        Ci[j*Cstride] = Ctmpi[j];
      }
    nfree(Ctmp);
}

void
simdotp(
  int M,
  int N,
  double *X,
  double *C,
  int weighted,
  double *W)
{
  time_t t0 = time(NULL);
  clock_t clock0 = clock();
  mess("computing pairwise correlation....\n");
  
  if( weighted == 0 )
    {
    #pragma omp parallel for schedule(runtime)
    for(int i = 0; i < N; i++ )
      {
      double *xi = X + i*M;
      double v = 0;
      for(int j = 0; j < M; j++ )
        v += xi[j]*xi[j];
      v = sqrt(v);
      if( v > 0 )
        for(int j = 0; j < M; j++ )
          xi[j] /= v;
      }

    xTx(M,N,X,1,C);
    }
  else if( weighted == 1 ) // weight = 0 if isnan() and 1 otherwise
    {
    // standardization
    #pragma omp parallel for schedule(runtime)
    for(int i = 0; i < N; i++ )
      {
      double *xi = X + i*M;
      double v = 0, w = 0;
      for(int j = 0; j < M; j++ )
        if( !isnan(xi[j]) )
          { v += xi[j]*xi[j]; w++; }
      v = sqrt(v);
      if( v > 0 )
        for(int j = 0; j < M; j++ )
          if( !isnan(xi[j]) ) xi[j] /= v;
      }

    double *WX = (double*)nalloc((size_t)M*N,sizeof(double));
    
    for(int ij = 0; ij < M*N; ij++ )
      WX[ij] = isnan(X[ij]) ? 0: X[ij];
    xTx( M, N, WX, 2, C );

    for(int ij = 0; ij < M*N; ij++ )
      WX[ij] = isnan(X[ij]) ? 0 : 1;
    xTx( M, N, WX, 2, C+1 );

    nfree(WX); 
    }
  else if( weighted == 2 )
    {
    // TODO: arbitrary weight matrix W
    }

  int elapsed = difftime( time(NULL), t0 );
  mess("elapsed : %4dh %2dm %2ds\n",
    elapsed/3600, (elapsed % 3600)/60, elapsed % 60);
  clock_t clock_now = clock();
  int cpu_t = round((clock_now - clock0)/CLOCKS_PER_SEC);
  mess("CPU time: %4dh %2dm %2ds\n",
    cpu_t/3600, (cpu_t % 3600)/60, cpu_t % 60);
}

void
R_sim(
  int* M_,
  int* N_,
  double *X,
  double *C,
  int* weighted_,
  double *W
)
{
  simdotp(*M_,*N_,X,C,*weighted_,W);
}
