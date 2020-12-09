/*
 * nnc.c - Nearest-Neighbor Chain Algorithm for Hierarchical Clustering
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "etc.h"
#include "nnc.h"

typedef struct
  {
  int id;
  double S;
  }
nnpq;

static void
nnpq_heapify( nnpq *q, int n )
{
  for(int i = 1; i < n; i++ )
    {
    nnpq tmp;
    if( q[i].S > q[0].S )
      {
      tmp = q[i];
      q[i] = q[0];
      q[0] = tmp;
      }
    
    tmp = q[i];
    int j = i;
    for(int k = j/2; k > 0 && (q[k].S > tmp.S); k /= 2 )
      {
      q[j] = q[k];
      j = k;
      }
    q[j] = tmp;
    }

  return;
}

static void
nnpq_push ( nnpq *q, int n, nnpq qin )
{
  if( n == 1 )
    {
    if( qin.S > q[0].S )
      q[0] = qin;
    return;
    }

  if( qin.S < q[1].S ) // ignore if smaller than the smallest
    return;  

  if( qin.S > q[0].S ) // swap if better than maximum
    {
    q[1] = q[0];
    q[0] = qin;
    qin = q[1];
    }

  int j = 1;
  for(int k = 2; k < n; k *= 2 )
    {
    if( (k+1 < n) && (q[k].S > q[k+1].S) )
      k++;

    if( qin.S <= q[k].S )
      break;

    q[j] = q[k];
    j = k;
    }
  q[j] = qin;

  return;
}

typedef struct
  {
  int id;
  int ncache;
  int next_merge;
  nnpq *q;
  }
nnchain;

int    // 0: ok, -1: memory allocation failure
nnc (

  // input data
  const int N,        // number of items
  void *data,         // data object, to be passed on to callback functions

  // callback functions
  void (*nnc_scan)(   // callback to compute scores between 
    void *data,
    int j,            // the anchor node
    int nK,           // number of nodes to compare
    int *K,           // list of node ids to compare to j
    double *Sjk       // list of S_{jk}
    ),

  void (*nnc_merge)(  // callback to merge node1 and node2 into new_node
    void *data,
    int new_node,
    int node1,
    int node2
    ),

  // tree output, storage allocated by caller
  int *L,             // [2*N] left child
  int *R,             // [2*N] right child
  int *U,             // [2*N] up (parent)
  double *S,          // [2*N] similarity score (branch height)
  
  // nnc specific options
  int cache_length,     // maximum cache length
  int verbose
  )
{
  /*
   *  allocate temporary workspace
   */
  nnchain *C = (nnchain*)nalloc( N, sizeof(nnchain) );
  if( NULL == C ) return -1;
  
  C[0].q = (nnpq*)nalloc( N*cache_length, sizeof(nnpq) );
  for(int i = 1; i < N; i++ )
    C[i].q = C[0].q + i*cache_length;
  
  int *tmp_id = (int*)nalloc(N,sizeof(int));
  double *tmp_S = (double*)nalloc(N,sizeof(double));

  if( NULL == C[0].q || NULL == tmp_id || NULL == tmp_S )
    { nfree(C[0].q); nfree(C); nfree(tmp_id); nfree(tmp_S); return -1; }

  /*
   *  initialize the free list 
   */
  for(int i = 1; i <= N; i++ ) // leaf nodes are in free list 
    {
    U[i] = 0;
    L[i] = i-1;
    R[i] = i+1;
    }
  // dummy at both ends of free list
  L[1] = 0; R[0] = 1;
  L[0] = N; R[N] = 0;

  // internal nodes
  for(int i = N+1; i < 2*N; i++ )
    U[i] = L[i] = R[i] = 0;

  /* 
  * The Algorithm
  */

  int freq_cache_emptied = 0;
  int freq_cache_hit = 0;
  int freq_scans = 0; // how often new whole-scans are performed
  int count_cmp = 0;  // count of pairwise score

#define PROGRESS_GRAIN 0.005 
  double progress = PROGRESS_GRAIN;

  time_t t0 = time(NULL);
  clock_t clock0 = clock();
  if( verbose )
    {
    mess("Number of items: %d\n\n",N );
    mess("%-15s%-20s%-10s%20s%14s\n",
        "time elapsed:","0%","   50%","100%","time left:");
    }

  int z = -1;  // end-of-chain pointer (empty)
  for(int next_merge = N+1; next_merge < 2*N; next_merge++)
    {
    // chain-extension loop
    while( R[0] != 0 ) // free list is not empty
      {
      if( z >= 0 )
        {
        if( C[z].ncache > 0 ) // refresh existing cache
          {
          // compact the nodes that are still not merged 
          int n = 0;
          for(int i = 0; i < C[z].ncache; i++ )
            if( U[ C[z].q[i].id ] == 0 )
              C[z].q[n++] = C[z].q[i];

          C[z].ncache = n;  // this can be zero!

          if( C[z].ncache > 0 ) 
            {
            nnpq_heapify( C[z].q, C[z].ncache );

            int n = 0;
            for(int i = C[z].next_merge; i < next_merge; i++ )
              if( U[i] == 0 )
                tmp_id[n++] = i;

            nnc_scan( data, C[z].id, n, tmp_id, tmp_S );
            count_cmp += n;

            for(int i = 0; i < n; i++ )
              {
              nnpq qi = { tmp_id[i], tmp_S[i] };
              nnpq_push( C[z].q, C[z].ncache, qi );
              }

            C[z].next_merge = next_merge;
            }
          else
            freq_cache_emptied++;
          }

        if( C[z].ncache == 0 ) // re-scan and rebuild the cache
          {
          int n = 0;
          for(int i = R[0]; i; i = R[i] )
            tmp_id[n++] = i;

          nnc_scan( data, C[z].id, n, tmp_id, tmp_S );
          freq_scans++;
          count_cmp += n;

          int p = n > cache_length ? cache_length : n;
          for(int i = 0; i < p; i++ )
            {
            C[z].q[i].id = tmp_id[i];
            C[z].q[i].S = tmp_S[i];
            }
          nnpq_heapify( C[z].q, p );
          C[z].ncache = p;

          for(int i = cache_length; i < n; i++ ) // skipped if n <= cache_length
            {
            nnpq qi = { tmp_id[i], tmp_S[i] };
            nnpq_push( C[z].q, cache_length, qi );
            }

          C[z].next_merge = next_merge; 
          }
        else
          freq_cache_hit++;
        }

      int new_tail = 0;

      if( z == -1 )
        new_tail = R[0];            // first item on free list
      else if( z == 0 )
        new_tail = C[z].q[0].id;    // NN to the chain end
      else  // z >= 1
        {
        if( C[z].q[0].S <= C[z-1].q[0].S ) // non-increasing link
          break;
        else
          new_tail = C[z].q[0].id;  // NN to the chain end
        }

      L[R[new_tail]] = L[new_tail];  // unlink from F
      R[L[new_tail]] = R[new_tail];
      L[new_tail] = R[new_tail] = 0;

      z++;
      C[z].id = new_tail;
      C[z].ncache = 0;

      } // chain-extension loop
    

    // merge the last two nodes on the chain
    int a = C[z-1].id, b = C[z].id;

    nnc_merge( data, next_merge, a, b );
    U[a] = U[b] = next_merge;
    S[next_merge] = C[z-1].q[0].S;
    
    // insert new node in front of the free node list
    R[next_merge] = R[0];
    L[R[next_merge]] = next_merge;
    R[0] = next_merge;
    L[next_merge] = 0;
    
    z -= 2; // pop

    double m = next_merge-N;
    double completed_frac = (m*(N-0.5*(m+1)))/((double)N*(N-1)/2);
    if( completed_frac > progress )
      {
      progress += PROGRESS_GRAIN;
      check_interrupt();
      if( verbose )
        {
        int elapsed = difftime( time(NULL), t0);
        mess("\r%4dh %2dm %2ds  ",
           elapsed/3600,(elapsed % 3600)/60,elapsed % 60 );
        for(int i = 0; i < 50; i++ )
          mess( "%c", i < floor(completed_frac*50) ? '#':'-');
        int remaining = round((elapsed/completed_frac)*(1-completed_frac));
        mess(" %4dh %2dm %2ds", 
            remaining / 3600, (remaining % 3600)/60, remaining % 60 );
        }
      }

    } // main loop

  // after the construction:
  if(verbose)
    {
    int elapsed = difftime( time(NULL), t0);
    mess("\r%4dh %2dm %2ds  ",
       elapsed/3600,(elapsed % 3600)/60,elapsed % 60 );
    mess("##################################################              \n");
    clock_t clock_now = clock();
    int cput = round((clock_now - clock0)/CLOCKS_PER_SEC);
    mess("%14s\n","CPU time:");
    mess("%4dh %2dm %2ds\n\n",
        cput/3600, (cput % 3600)/60, cput % 60 );
    mess("cache_length = %d\nrescans      = %d\ncache hits   = %d\n"
         "cache empty  = %d\ncomparisons  = %d\n\n",
        cache_length,
        freq_scans, freq_cache_hit, freq_cache_emptied, count_cmp );
    }

  /*
   *  fill-in links to child nodes
   */
  for(int i = 1; i < 2*N-1; i++ )
    {
    int j = U[i];
    if( L[j] == 0 )  // left child has smaller id
      L[j] = i;
    else
      R[j] = i;
    }

  /*
   * free workspace
   */
  nfree(tmp_S);
  nfree(tmp_id);
  nfree(C[0].q);
  nfree(C);
  return 0;
}

