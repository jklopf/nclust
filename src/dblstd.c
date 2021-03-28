/*
 * dblstd.c
 *
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void
dblstd (
  const int *dims,
  double *y,
  int *iter_max_,
  double *tol_,
  int *verbose_
  )
{
  const int m = dims[0], n = dims[1];
  double *s0 = (double*) malloc(sizeof(double) * m * 2 );
  double *s1 = s0 + m;

  for(int i = 0; i < m; i++ )
    s0[i] = 1;

  int iter;

  for(iter = 1; iter <= iter_max_[0]; iter++ )
    {
    for(int i = 0; i < m; i++ )
      s1[i] = 0;

    for(int j = 0; j < n; j++ )
      {
      double *yj = y + j * m;
      double sj = 0;
      #pragma omp simd
      for(int i = 0; i < m; i++ )
        {
        yj[i] /= s0[i];
        sj += yj[i]*yj[i];
        }
      sj = sqrt( sj/m );
      if( sj > 0 )
        for(int i = 0; i < m; i++ )
          {
          yj[i] /= sj;
          s1[i] += yj[i]*yj[i];
          }
      }

    double S = 0;
    for(int i = 0; i < m; i++ )
      {
      s0[i] = sqrt( s1[i]/n );
      if(s1[i] == 0 ) s0[i] = 1;
      double ds = fabs(s0[i]-1);
      if( ds > S ) S = ds;
      }
    if(*verbose_) 
      fprintf(stderr,"dblstd: %3d %g\n",iter,S);
    if( S < *tol_ )
      {
      *tol_ = S;
      break;
      }
    }

  iter_max_[0] = iter;
  free(s0);
}
