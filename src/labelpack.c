#include <stdlib.h>
#include <stdio.h>

#include "etc.h"

// arrange placement of sparse labels to minimize overlap
//
// TODO: this should be a part of generic graphic utilities
//
void
R_labelpack (
  int *n_, const double *x,   // original positions, assume sorted already
  double *d_,                  // space taken by each, same unit as x
  double *minx_, double *maxx_, // boundary constraints on the output
  int *maxiter_,                // maximum iteration
  double *y                     // output
  )
{
  int n = *n_;
  if(n <= 0) return;

  double d = *d_;
  double minx = *minx_;
  double maxx = *maxx_;

  for(int i = 0; i < n; i++ )
    {
    y[i] = x[i];
    if(y[i] < minx ) y[i] = minx;
    if(y[i] > maxx ) y[i] = maxx;
    }
    
  int *g = (int*)nalloc(n+1, sizeof(int));
  for(int i = 0; i <= n; i++ )
    g[i] = i;

  int collision = 0;
  for(int r = 0; r < *maxiter_; r++ )
    {
    collision = 0;
    for( int i = g[n-1]; i > 0; )
      {
      if( y[i] - y[i-1] <= d )
        {
        collision = 1;
        for(int j = i; g[j] == i; j++ )
          g[j] = g[i-1];
        }
      i = g[i-1];
      }
      
    if(!collision) break;
    
    /* adjust the location */
    for(int i = n-1; i > 0; i = g[i]-1 )
      {
      int nrun = i-g[i]+1;
      if(nrun == 1) continue;
      double mean_x = 0;
      for(int j = g[i]; j <= i; j++ )
        mean_x += y[j];
      mean_x /= nrun;
      double mean_j = 0.5*(g[i] + i);
      for(int j = g[i]; j <= i; j++ )
        y[j] = mean_x + (j-mean_j)*d;
        
      if(y[g[i]] < minx && y[i] > maxx )
        mess("overcrowding\n");
        
      else if(y[g[i]] < minx )
        {
        double u = minx - y[g[i]] ;
        for(int j = g[i]; j <= i; j++ )
          y[j] += u;
        }
      else if( y[i] > maxx )
        {
        double u = y[i] - maxx;
        for(int j = g[i]; j <= i; j++ )
          y[j] -= u;
        }
      }
    }

  nfree(g);
  return;
}

