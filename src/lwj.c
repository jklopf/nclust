/*
 * lwj.c 
 */

/*
INTERFACE:

INPUT: 
- a half-matrix, N*(N+1)/2 contiguous in memory of pairwise
  score between leaf nodes (to be overwritten by the workspace)
- linkage: (default: average)
- option switches: cache, verbose

OUTPUT:
- tree structure (as defined in `nnc.h`)


IMPLEMENTATION:
CALLBACKS:

- scan: 
  scan the half-matrix (no need to compute!), but need to
  do indexing work. Note that nnc define its node space
  in {0,1,...,N,N+1,...,2*N-1}. lwj needs to define
  the mapping to re-cycled matrix index.

- merge:
  Lance-William recurrence to combine two rows/columns into
  one. The two original rows' values are not needed anymore
  and thus the lower index can be overwritten, while the
  higher one will be effectively ignored after the merge.
  The three node indices need to be converted to matrix indices.
  The new node is index to a table that lwj keep, to map to
  recycled matrix index.

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nnc.h"
#include "etc.h"

typedef struct
  {
  int N;
  double **C_;  // [N] pointers to the base of column vectors
  int *imap;    // [2*N] map node index to matrix index
  double *w;    // [N] node weights
  int nact;     // number of active matrix index
  int *acti;    // [N] list of active matrix index, terminated by -1
  int *iact;    // [N] convert matrix index to `active` list position
  int method;   // link method
  int weighted; // if C_ has weights
  }
lwj_context;

void lwj_scan(
  void *ctx_,
  int j,
  int nK,
  int *K,
  double *Sjk
  )
{
  lwj_context *ctx = (lwj_context*)ctx_;
  int *imap = ctx->imap;
  double **C_ = ctx->C_;
  int weighted = ctx->weighted;
  static int counter=0;

  int a = imap[j];
  if( !weighted )
    for(int i = 0; i < nK; i++ )
      {
      int b = imap[K[i]];
      Sjk[i] = ( a < b ? C_[b][a] : C_[a][b] );
      }
  else
    for(int i = 0; i < nK; i++ )
      {
      int b = imap[K[i]];
      Sjk[i] = ( a < b ? C_[b][2*a] : C_[a][2*b] );
      }
}

void lwj_merge(
  void *ctx_,
  int new_node,
  int node1,
  int node2
  )
{
  lwj_context *ctx = (lwj_context*)ctx_;
  int N = ctx->N;
  int *imap = ctx->imap;
  double **C_ = ctx->C_;
  double *w = ctx->w;
  int *acti = ctx->acti;
  int *iact = ctx->iact;
  int method = ctx->method;

  int a, b;
  if( imap[node1] < imap[node2] )
    { a = imap[node1]; b = imap[node2]; }
  else
    { a = imap[node2]; b = imap[node1]; }

  if( ctx->weighted == 0 )
    {
    double Sw_ab = w[a] + w[b];
    switch(method)
      {
      case 0: // average, aka UPGMA
      default:
        #pragma omp parallel for schedule(runtime)
        for(int i = 0; i < iact[a]; i++ )
          {
          int c = acti[i];
          C_[a][c] = (w[a]*C_[a][c] + w[b]*C_[b][c])/Sw_ab; 
          }
        #pragma omp parallel for schedule(runtime)
        for(int i = iact[a]+1; i < iact[b]; i++ )
          {
          int c = acti[i];
          C_[c][a] = (w[a]*C_[c][a] + w[b]*C_[b][c])/Sw_ab;
          }
        #pragma omp parallel for schedule(runtime)
        for(int i = iact[b]+1; i < ctx->nact; i++ )
          {
          int c = acti[i];
          C_[c][a] = (w[a]*C_[c][a] + w[b]*C_[c][b])/Sw_ab;
          }
        break;

      case 1: // MISSQ (a.k.a Ward; same as`ward.D` in R, not `ward.D2`)
        #pragma omp parallel for schedule(runtime)
        for(int i = 0; i < iact[a]; i++ )
          {
          int c = acti[i];
          double d = w[c]*C_[b][a];
          double Sw = Sw_ab + w[c];
          C_[a][c] = ((w[a]+w[c])*C_[a][c] + (w[b]+w[c])*C_[b][c] - d)/Sw;
          }
        #pragma omp parallel for schedule(runtime)
        for(int i = iact[a]+1; i < iact[b]; i++ )
          {
          int c = acti[i];
          double d = w[c]*C_[b][a];
          double Sw = Sw_ab + w[c];
          C_[c][a] = ((w[a]+w[c])*C_[c][a] + (w[b]+w[c])*C_[b][c] - d)/Sw;
          }
        #pragma omp parallel for schedule(runtime)
        for(int i = iact[b]+1; i < ctx->nact; i++ )
          {
          int c = acti[i];
          double d = w[c]*C_[b][a];
          double Sw = Sw_ab + w[c];
          C_[c][a] = ((w[a]+w[c])*C_[c][a] + (w[b]+w[c])*C_[c][b] - d)/Sw;
          }
        break;

      case 2: // single 
        #pragma omp parallel for schedule(runtime)
        for(int i = 0; i < iact[a]; i++ )
          {
          int c = acti[i];
          C_[a][c] = (C_[a][c] + C_[b][c] + fabs(C_[a][c] - C_[b][c]) )/2; 
          }
        #pragma omp parallel for schedule(runtime)
        for(int i = iact[a]+1; i < iact[b]; i++ )
          {
          int c = acti[i];
          C_[c][a] = (C_[c][a] + C_[b][c] + fabs(C_[c][a] - C_[b][c]))/2;
          }
        #pragma omp parallel for schedule(runtime)
        for(int i = iact[b]+1; i < ctx->nact; i++ )
          {
          int c = acti[i];
          C_[c][a] = (C_[c][a] + C_[c][b] + fabs(C_[c][a] - C_[c][b]))/2;
          }
        break;

      case 3: // complete
        #pragma omp parallel for schedule(runtime)
        for(int i = 0; i < iact[a]; i++ )
          {
          int c = acti[i];
          C_[a][c] = (C_[a][c] + C_[b][c] - fabs(C_[a][c] - C_[b][c]) )/2; 
          }
        #pragma omp parallel for schedule(runtime)
        for(int i = iact[acti[iact[a]+1]]; i < iact[b]; i++ )
          {
          int c = acti[i];
          C_[c][a] = (C_[c][a] + C_[b][c] - fabs(C_[c][a] - C_[b][c]))/2;
          }
        #pragma omp parallel for schedule(runtime)
        for(int i = iact[acti[iact[b]+1]]; i < ctx->nact; i++ )
          {
          int c = acti[i];
          C_[c][a] = (C_[c][a] + C_[c][b] - fabs(C_[c][a] - C_[c][b]))/2;
          }
        break;
      }
    w[a] = Sw_ab; 
    }
  else // weighted: only ave link is sensible
    {
    int ac = 2*a, aw = 2*a+1;
    int bc = 2*b, bw = 2*b+1;

    #pragma omp parallel for schedule(runtime)
    for(int i = 0; i < iact[a]; i++ )
      {
      int c = acti[i];
      int cc = 2*c, cw = 2*c+1;
      double Sw = C_[a][cw] + C_[b][cw];
      if( Sw == 0 ) 
        { C_[c][ac] = C_[c][aw] = 0; continue; }
      C_[a][cc] = (C_[a][cw]*C_[a][cc] + C_[b][cw]*C_[b][cc])/Sw;
      C_[a][cw] = Sw;
      }
    #pragma omp parallel for schedule(runtime)
    for(int i = iact[a]+1; i < iact[b]; i++ )
      {
      int c = acti[i];
      int cc = 2*c, cw = 2*c+1;
      double Sw = C_[c][aw] + C_[b][cw]; 
      if( Sw == 0 ) 
        { C_[c][ac] = C_[c][aw] = 0; continue; }
      C_[c][ac] = (C_[c][aw]*C_[c][ac] + C_[b][cw]*C_[b][cc])/Sw;
      C_[c][aw] = Sw;
      }
    #pragma omp parallel for schedule(runtime)
    for(int i = iact[b]+1; i < ctx->nact; i++ )
      {
      int c = acti[i];
      int cc = 2*c, cw = 2*c+1;
      double Sw = C_[c][aw] + C_[c][bw];
      if( Sw == 0 ) 
        { C_[c][ac] = C_[c][aw] = 0; continue; }
      C_[c][ac] = (C_[c][aw]*C_[c][ac] + C_[c][bw]*C_[c][bc])/Sw;
      C_[c][aw] = Sw;
      }
    }

  (ctx->nact)--;
  for(int i = iact[b]; i < ctx->nact; i++ )
    {
    acti[i] = acti[i+1];
    iact[acti[i]] = i;
    }
  iact[b] = -1;
  
  imap[new_node] = a;
}

int
lwj (
  int N,     // number of leaf nodes
  double *C, // [N*(N+1)/2 * (weighted?1:2)] similarity matrix, diagonal unused

  int *L,     // [2*N] left
  int *R,     // [2*N] right
  int *U,     // [2*N] up (parent)
  double *S,  // [2*N] branch height

  int method, // 0: average, 1: Ward
  int weighted, // 0: no, !=0: yes
  int max_ncache,
  int verbose
  )
{
  // init context
  lwj_context ctx;

  ctx.N = N;

  // pointers to beginning of column vectors
  ctx.C_ = (double**) nalloc(N,sizeof(double*));
  for(size_t i = 0; i < N; i++ )
    ctx.C_[i] = C + (i*(i+1)/2) * (weighted?2:1);

  // map node index to matrix index
  ctx.imap = (int*) nalloc( 2*N, sizeof(int) );
  ctx.imap[0] = 0;
  for(int i = 1; i <= N; i++ )
    ctx.imap[i] = i-1;
  for(int i = N+1; i < 2*N; i++ )
    ctx.imap[i] = 0;

  // active matrix index list
  ctx.nact = N;
  ctx.acti = (int*)nalloc(N,sizeof(int) ); // map matrix id 0:(N-1) to list id
  ctx.iact = (int*)nalloc(N,sizeof(int) ); // reverse of the above
  for(int i = 0; i < N; i++ )
    {
    ctx.acti[i] = i;
    ctx.iact[i] = i;
    }

  ctx.method = method;
  ctx.weighted = weighted;

  // marginal node weights (only use in unweighted problems)
  ctx.w = (double*) nalloc( N, sizeof(double) );
  for(int i = 0; i < N; i++ )
    ctx.w[i] = 1;

  nnc( N, &ctx,
    &lwj_scan,
    &lwj_merge,
    L,R,U,S,
    max_ncache, verbose );

  nfree(ctx.C_);
  nfree(ctx.imap);
  nfree(ctx.acti);
  nfree(ctx.iact);
  nfree(ctx.w);
}
