/*
 * wdownsample.c - subset a matrix, tile, and compute the mean of each tile 
 *
 * - subset is specified by integer indices, which may contain
 *   replicate selections, or missing values (represented by -1)
 *   and may be arbitrarily re-ordered
 * - averaging is done over integer-sized windows for rows and columns
 * - the last window is padded by missing indices
 * - caller responsible in ensuring irow and icol are within 0:(mx-1)
 *   and 0:(nx-1)
 */

#include <stdio.h>

#include <math.h>

#define NULL_ID -1   // index to nothing

void
wdownsample (
  int *dims, // mx, nx, my, ny, mw, nw,
  int *irow, // row indices (of length my * mw)
  int *icol, // col indices (of length ny * nw)
  double *x, // input matrix of size mx * nx [column vector is contiguous]
  double *wx, // input weights
  double *y,  // output matrix of size my * ny
  double *wy, // output weights
  double *wrow,  // sum of non-NA indices for each tile row
  double *wcol   // sum of non-NA indices for each tile column
  )
{
  int mx = dims[0];  // input size
  int nx = dims[1];
  int my = dims[2];  // output size
  int ny = dims[3];
  int mw = dims[4];  // window size
  int nw = dims[5];

  for(int i = 0; i < my; i++ )
    {
    wrow[i] = 0;
    for(int u = 0; u < mw; u++ )
      if( irow[i*mw + u] != NULL_ID )
        wrow[i]++;
    }

  for(int j = 0; j < ny; j++ )
    {
    wcol[j] = 0;
    for(int v = 0; v < nw; v++ )
      if( icol[j*nw + v] != NULL_ID  )
        wcol[j]++;
    }

  for(int j = 0; j < ny; j++ )
    {
    if( wcol[j] == 0 )
      {
      for(int i = 0; i < my; i++ )
        {
        y[i + j*my] = 0;
        wy[i + j*my] = 0;
        }
      continue;
      }
    for(int i = 0; i < my; i++ )
      {
      if( wrow[i] == 0 )
        {
        y[i + j*my] = 0;
        wy[i + j*my] = 0;
        continue;  
        }

      double sy = 0;
      double sw = 0;

      for(int v = 0; v < nw; v++ )
        {
        int jx = icol[j*nw + v];
        if( jx == NULL_ID ) continue;
        for(int u = 0; u < mw; u++ )
          {
          int ix = irow[i*mw + u];
          if(ix == NULL_ID ) continue;
          double w = wx[ ix + jx * mx ];
          sy += w * x[ ix + jx * mx ];
          sw += w;
          }
        }
      
      wy[i + j*my] = sw;
      y[i + j*my] = sw > 0 ?  sy/sw: 0;
      }
    }
  return;
}
