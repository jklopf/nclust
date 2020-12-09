#include <stdlib.h>
#include <math.h>

#include "sim.h"
#include "lwj.h"
#include "nclustree.h"

#include "etc.h"

void
R_lwj(
  int *M_,
  int *N_,
  int *weighted_,  // if 0 not weighted
  double *X,    // [MN]
  double *W,    // [MN, optional]
  int *L,
  int *R,
  int *U,
  double *S,
  int *order,
  int *nleaf,
  int *leftmost,
  int *level,
  const int *method_,
  const int *flip_,
  const int *verbose_
  )
{
  int N = *N_, M = *M_;
  int weighted = *weighted_;
  int method = *method_;
  
  double *C = (double*)nalloc( (size_t)N * (N+1)/2 * (weighted ? 2 : 1),
    sizeof(double) );

  simdotp(M,N,X,C,weighted,W);

  mess("\ncomputing hierarchical clustering...\n");
  lwj( N, C, L, R, U, S, method, weighted, 32, *verbose_ );
  nfree(C);

  /********************* post clustering ******************************/
  // TODO: put into separate function so we don't
  // have to redo the cov/cor matrix and lwj.
  mess("branch flip = %d\n",*flip_);
  switch(*flip_)
    {
    case 0:
    default:
      break;
    case 1:
      branchflip_tightleft( N, U, L, R, S );
      break;
    case 2:
      branchflip_tightright( N, U, L, R, S );
      break;
    case 3:
      if( weighted == 0 )
        {
        double *Y = (double*)nalloc((size_t)M*N,sizeof(double));
        double *wnode = (double*)nalloc((size_t)2*N,sizeof(double));
        for(int i = 1; i <= N; i++ )
          wnode[i] = 1.0;
        for(int i = N+1; i < 2*N; i++ )
          wnode[i] = 0.0;
        for(int i = 1; i < 2*N-1; i++ )
          {
          double *Yi = (i <= N ? X + (i-1)*M : Y + (i-N-1)*M);
          int u = U[i];
          double *Yu = Y + (u-N-1)*M;
          for(int k = 0; k < M; k++ )
            Yu[k] = (wnode[u]*Yu[k] + wnode[i]*Yi[k])/(wnode[u]+wnode[i]);  
          wnode[u] += wnode[i];
          }

        double bflip_func( void *data, int i, int j)
          {
          double *Yi = (i <= N ? X + (i-1)*M : Y + (i-N-1)*M );
          double *Yj = (j <= N ? X + (j-1)*M : Y + (j-N-1)*M );
          double v = 0;
          for(int k = 0; k < M; k++ )
            v += Yi[k] * Yj[k];
          return v;
          }

        branchflip_nnephew( U, L, R, NULL, bflip_func );

        nfree(wnode);
        nfree(Y);
        }
      else if(weighted == 1 )
        {
        double *Y = (double*)nalloc((size_t)M*N*2,sizeof(double));
        double *wY = (double*)nalloc((size_t)M*N*2,sizeof(double));
        for(int ij = 0; ij < M*N; ij++ )
          {
          Y[M+ij] = isnan(X[ij]) ? 0: X[ij];
          wY[M+ij] = isnan(X[ij]) ? 0: 1;
          }
        for(int ij = M*(N+1); ij < M*N*2; ij++ )
          wY[ij] = Y[ij] = 0;
        for(int i = 1; i < 2*N; i++ )
          {
          double *Yi = Y + i*M, *Yu = Y + U[i]*M;
          double *wYi = wY + i*M, *wYu = wY + U[i]*M;
          for(int k = 0; k < M; k++ )
            {
            double sw = wYu[k] + wYi[k];
            if( sw > 0 )
              Yu[k] = (wYu[k]*Yu[k] + wYi[k]*Yi[k])/sw;
            wYu[k] = sw;
            }
          }
        double bflip_func2( void *data, int i, int j )
          {
          double *Yi = Y + i*M, *Yj = Y + j*M;
          double *wYi = wY + i*M, *wYj = wY + j*M;
          double v = 0, w = 0;
          for(int k = 0; k < M; k++ )
            {
            w += wYi[k] * wYj[k];
            v += wYi[k]*Yi[k] * wYj[k] * Yj[k];
            }
          //mess("%g %g %g\n",v, w,v/w);
          return w > 0 ? v/w : 0;
          }
        branchflip_nnephew( U, L, R, NULL, bflip_func2 );

        nfree(Y); nfree(wY);
        }
      break;
    }
  branch_nleaf( U, L, R, nleaf );
  leafordering( U, L, R, order, leftmost );
  branch_level( &N, U, level );
}
