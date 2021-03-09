/*
 * wf_dense.c - dense weighted matrix data structure for nnc
 *
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "nnc.h"
#include "nclustree.h"
#include "wf_dense.h"

//
// standardize 
//
static void
wf_dense_stdz (
  wf_dense *D
  )
{
  if( D->w )
    {
    #pragma omp parallel for schedule(runtime)
    for(int j = 1; j <= D->N; j++ )
      {
      double *xj = D->wx + D->M * j;
      double *wj = D->w + D->M * j;
      double sw = 0, swx2 = 0;
      for(int i = 0; i < D->M; i++ )
        {
        double ww = D->t[i]*D->t[i]*wj[i]*wj[i];
        sw += ww;
        swx2 += ww * xj[i] * xj[i];
        }
      if( swx2 > 0 )
        {
        double v = swx2 / sw;
        for(int i = 0; i < D->M; i++ )
          wj[i] *= v;
        double sd = sqrt( v );
        for(int i = 0; i < D->M; i++ )
          xj[i] /= sd;
        }
      }
    }
  else
    {
    #pragma omp parallel for schedule(runtime)
    for(int j = 1; j <= D->N; j++ )
      {
      double *xj = D->wx + D->M * j;
      double sw = D->St2, swx2 = 0;
      for(int i = 0; i < D->M; i++ )
        {
        swx2 += D->t[i] * D->t[i] * xj[i] * xj[i];
        }
      if( swx2 > 0 )
        {
        double sd = sqrt( swx2 / sw );
        for(int i = 0; i < D->M; i++ )
          xj[i] /= sd;
        }
      }
    }
}

// if w is null, wx contains t_i * x
// if w is not null, wx = t_i * w_ij * x
static void
wf_dense_premultiply_avedot (
  wf_dense *D
  )
{
  if( D->w )
    {
    #pragma omp parallel for schedule(runtime)
    for(int j = 1; j <= D->N; j++ )
      {
      double *wxj = D->wx + D->M * j;
      double *wj = D->w + D->M * j;
      for(int i = 0; i < D->M; i++ )
        {
        wj[i] *= D->t[i];
        wxj[i] *= wj[i];
        }
      }
    }
  else
    {
    #pragma omp parallel for schedule(runtime)
    for(int j = 1; j <= D->N; j++ )
      {
      double *wxj = D->wx + D->M * j;
      for(int i = 0; i < D->M; i++ )
        wxj[i] *= D->t[i];
      }
    }
}

static void
wf_dense_premultiply_ward (
  wf_dense *D
  )
{
  if( D->w )
    {
    #pragma omp parallel for schedule(runtime)
    for(int j = 1; j <= D->N; j++ )
      {
      double *wxj = D->wx + D->M * j;
      double *wj = D->w + D->M * j;
      for(int i = 0; i < D->M; i++ )
        wj[i] *= D->t[i] * D->u[j];
      }
    }
}

// 
// find the nearest for a given node j
//
void
wf_dense_nnc_scan_avedot(
  void *data,
  int j,
  int nc,
  int *c_,
  double *Sj_
  )
{
  wf_dense *D = (wf_dense*)data;
  const int M = D->M;

  const double *wxj = D->wx + M*j;
  #pragma omp parallel for schedule(runtime)
  for(int k = 0; k < nc; k++ )
    {
    const double *wxk = D->wx + M*c_[k];
    double swx = 0;
    #pragma omp simd
    for(int i = 0; i < M; i++ )
      swx += wxj[i] * wxk[i];
    Sj_[k] = swx;
    }

  if( D->w )
    {
    const double *wj = D->w + M*j;
    #pragma omp parallel for schedule(runtime)
    for(int k = 0; k < nc; k++ )
      {
      const double *wk = D->w + M*c_[k];
      double sw = 0;
      #pragma omp simd
      for(int i = 0; i < M; i++ )
        sw += wj[i] * wk[i];
      if(sw > 0 ) Sj_[k] /= sw;
      else Sj_[k] = 0; // zero is indifference for cov & cor
      }
    }
  else
    {
    for(int k = 0; k < nc; k++ )
      {
      Sj_[k] /= D->St2;
      }
    }

  return;
}

// 
// Ward's method
//
void
wf_dense_nnc_scan_ward(
  void *data,
  int j,
  int nc,
  int *c_,
  double *Sj_
  )
{
  wf_dense *D = (wf_dense*)data;
  const int M = D->M;

  const double *xj = D->wx + M*j;

  if( D->w )
    {
    const double *wj = D->w + M*j;

    #pragma omp parallel for schedule(runtime)
    for(int k = 0; k < nc; k++ )
      {
      const double *xk = D->wx + M*c_[k];
      const double *wk = D->w + M*c_[k];
      double sxx = 0;
      #pragma omp simd
      for(int i = 0; i < M; i++ )
        {
        double d = xj[i] - xk[i];
        double sumw = wj[i] + wk[i];
        double ww = (wj[i]*wk[i]);
        if( sumw > 0 ) ww /= sumw;
        sxx -= D->t[i]*ww*d*d;
        }
      Sj_[k] = sxx;
      }
    }
  else
    {
    double wj = D->u[j];

    #pragma omp parallel for schedule(runtime)
    for(int k = 0; k < nc; k++ )
      {
      const double *xk = D->wx + M*c_[k];
      double sxx = 0;
      #pragma omp simd
      for(int i = 0; i < M; i++ )
        {
        double d = xj[i] - xk[i];
        sxx -= D->t[i]*d*d;
        }
      double wk = D->u[c_[k]];
      Sj_[k] = wj*wk/(wj+wk) * sxx;
      }
    }

  return;
}

void
wf_dense_nnc_merge_avedot(
  void *data,
  int m,
  int j,
  int k
  )
{
  wf_dense *D = (wf_dense*)data;
  const int M = D->M;
  double uj = D->u[j];
  double uk = D->u[k];
  double um = uj + uk;
  D->u[m] = um;

  double *wxj = D->wx + M*j;
  double *wxk = D->wx + M*k;
  double *wxm = D->wx + M*m;

  for(int i = 0; i < M; i++ )
    wxm[i] = (uj * wxj[i] + uk * wxk[i])/um;

  if( D->w )
    {
    double *wj = D->w + M*j;
    double *wk = D->w + M*k;
    double *wm = D->w + M*m;
    
    for(int i = 0; i < M; i++ )
      wm[i] = (uj * wj[i] + uk * wk[i])/um;
    }
}

void
wf_dense_nnc_merge_ward (
  void *data,
  int m,
  int j,
  int k
  )
{
  wf_dense *D = (wf_dense*)data;
  const int M = D->M;
  double *xj = D->wx + M*j;
  double *xk = D->wx + M*k;
  double *xm = D->wx + M*m;

  if( D->w )
    {
    double *wj = D->w + M*j;
    double *wk = D->w + M*k;
    double *wm = D->w + M*m;

    for(int i = 0; i < M; i++ )
      {
      wm[i] = wj[i] + wk[i];
      if(wm[i] > 0)
        xm[i] = (wj[i]*xj[i] + wk[i]*xk[i])/wm[i];
      else
        xm[i] = 0;
      }
    }
  else
    {
    D->u[m] = D->u[j] + D->u[k];
    for(int i = 0; i < M; i++ )
      xm[i] = (D->u[j]*xj[i] + D->u[k]*xk[i])/D->u[m];
    }
}



void
wf_dense_nclust
  (
  const int *dims,
  const int *options,

  double *xx,   // feature data matrix + workspace for branches
  double *ww,   // joint weights + workspace for branches
  double *t,    // feature weights
  double *u,    // item weights

  int *L, int *R, int *U,
  double *S,
  int *order,
  int *nleaf,
  int *leftmost,
  int *level,
  int *retstat
  )
{
  const int W = dims[0];  // 1: ww is null, 2 otherwise
  const int M = dims[1];
  const int N = dims[2];
  const int cache_length = options[0];
  const int branchflip = options[1];
  const int standardize = options[2];
  const int verbose = options[3];
  const int method = options[4];

  if( W == 1 )
    ww = NULL;

  double St2 = 0;
  for(int i = 0; i < M; i++ )
    St2 += t[i]*t[i];

  wf_dense data = { M, N, xx, ww, t, u, St2 };

  // standardize the data
  if( standardize )
    wf_dense_stdz( &data);

  // pre-multiply the weights
  if( method == 1 )
    wf_dense_premultiply_avedot(&data);

  // 
  // The Thing
  //
  nnc( N, &data, 
    method == 1 ? wf_dense_nnc_scan_avedot : wf_dense_nnc_scan_ward,
    method == 1 ? wf_dense_nnc_merge_avedot : wf_dense_nnc_merge_ward,
    L, R, U, S, cache_length, verbose );
  
  // branch flipping method. Default to the original NN-chain order.
  if( branchflip == 1 ) // centered by nearest nephew
    {
    if( S[ L[R[0]] ] < S[ R[R[0]] ] ) // arbitrary rule for the root
      SWAP(int, L[R[0]], R[R[0]] );

    branchflip_nnephew( U, L, R, &data,
     method == 1 ?  wf_dense_wcov : wf_dense_ward );
    }
  else if( branchflip == 2 )
    branchflip_tightleft( N, U, L, R, S );
  else if( branchflip == 3 ) 
    branchflip_tightright( N, U, L, R, S );

  // tree auxiliary info
  // - just for convenience, as they can all be re-derived from U, L and R
  // 
  branch_nleaf( U, L, R, nleaf );
  leafordering( U, L, R, order, leftmost );
  branch_level( &N, U, level );

  *retstat = 0;
  return;
}
