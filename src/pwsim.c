/*
 * pwsim.c - pairwise similarity matrix
 *
 */

#include <math.h>
#include "pwsim.h"

// dense pairwise dot-product
void
pwdot(
  const int *M_,
  const int *N_,
  const double *X,
  double *C           // [N(N+1)/2] upper triangular
  )
{
  int M = *M_, N = *N_;
  const double *xi = X;
  double *Ci = C;
  for(int i = 0; i < N; i++, xi += M, Ci += i )
    {
    #pragma omp parallel for schedule(runtime)
    for(int j = 0; j <= i; j++ )
      {
      const double *xj = X + j*M;
      double v = 0;
      for(int k = 0; k < M; k++ )
        v += xi[k]*xj[k];
      Ci[j] = v;
      }
    }
}

// scale by the product of the sqrt(diag(C))
void
pwscale(
  const int *N_,
  double *C
  )
{
  const int N = *N_;
  double *Ci = C;
  for(int i = 0; i < N; i++, Ci += i )
    {
    Ci[i] = sqrt(Ci[i]);
    for(int j = 0; j < i; j++ )
      Ci[j] /= Ci[i];

    double *Cji = Ci + 2*i + 1;
    for(int ii = i+1; ii < N; Cji += ++(ii) )
      *Cji /= Ci[i];
    }
}
