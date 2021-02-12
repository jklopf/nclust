/*
 * nclustree.c
 *
 * Various functions for handling nclust tree data structure.
 *
 * The binary tree is represented by three integer index links per node:
 * U(p), L(eft) and R(ight). This is a doubly-linked binary tree.
 *
 * - The zero-eth node is a dummy node.
 * - L[0] == R[0] is the root. U[0] == 0.
 * - The root point to zero as the parent.
 * - The leaves point to zero as their left and right children.
 *
 * Empty and singleton tree is possible.
 *
 * For a tree with N nodes, the length of link arrays is 2*N.
 *
 * nclust uses index 1:N for leaf nodes, and the root is always the
 * last node (index: 2*N-1).
 *
 * Depth-first tree traversal can be done with neither recursion nor stack.
 * The doubly-linked tree allows state-machine
 * iteration using a "stack" implicit in the links.
 * The `state` keeps track the stage of the visits.
 *
 */
#include <stdlib.h>
#include <stdio.h> // debugging only

#include "nclustree.h"

#ifndef SHELLPROG 
  #include <R.h> // for Rprintf
  #define mess(...) Rprintf(__VA_ARGS__)
#else
  #define mess(...) fprintf(stderr,__VA_ARGS__)
#endif


/* Given a binary tree, return the ordering of the leaves and
 * the left-most leaf of each cluster.
 */
void
leafordering( 
  const int *U,
  const int *L,
  const int *R,
  int *order,
  int *leftmost )
{
  int r = 0;    // index to 'order'
  int p = 0; // previous node

  for(int i = R[0]; i; )
    {
    if( p == U[i] ) // arrive from parent
      {
      if( L[i] == 0 ) // in a leaf node
        {
        order[r] = i;
        leftmost[i] = r;
        r++;
        
        p = i;
        i = U[i]; // return to parent
        }
      else
        {
        p = i;
        i = L[i]; // visit left node
        }
      }
    else if( p == L[i] ) // coming back from the left
      {
      p = i;
      i = R[i];  // visit the right node
      }
    else // coming back from the right
      {
      leftmost[i] = leftmost[ L[i] ];
      
      p = i;
      i = U[i]; // return to parent
      }
    }
}

/* count the number of leaves under each node
 *
 * - leaf nodes has count equals 1
 */
void
branch_nleaf (
  const int *U,
  const int *L,
  const int *R,
  int *nleaf
  )
{
  nleaf[0] = 0;
  for(int i = 1; i <= L[0]; i++ )
    if(L[i] == 0 && R[i] == 0 )
      nleaf[i] = 1;
    else
      nleaf[i] = nleaf[L[i]]+nleaf[R[i]];
}

/*
 * The level of the branches, starting with 1 at the root
 */

void
branch_level (
  const int *N_,
  const int *U,
  int* level
  )
{
  int N = *N_;
  level[2*N-1] = 1;
  for(int i = 2*N-2; i; i-- )
    level[i] = level[U[i]] + 1;
}

void
branchflip_tightleft(
  int N,
  const int *U,
  int *L,
  int *R,
  const double *S
  )
{
  for(int i = N+1; i < 2*N; i++ )
    if( L[i] <= N && R[i] <= N && L[i] > R[i]
      || L[i] > N && R[i] > N && S[L[i]] < S[R[i]]
      || L[i] > N && R[i] <= N )
      SWAP(int,L[i],R[i]);
}

void
branchflip_tightright(
  int N,
  const int *U,
  int *L,
  int *R,
  const double *S
  )
{
  for(int i = N+1; i < 2*N; i++ )
    if( L[i] <= N && R[i] <= N && L[i] > R[i]
      || L[i] > N && R[i] > N && S[L[i]] > S[R[i]]
      || L[i] > N && R[i] <= N )
      SWAP(int,L[i],R[i]);
}

/*
 * nearest-nephew branch flipping method
 *
 * It's a top-down traversal, where the flip is done at the top first.
 * It assumes the similarity between internal nodes can be calculated
 * efficiently at the top.
 */
void
branchflip_nnephew(
  const int *U,
  int *L,
  int *R,

  void *data,
  double (sim)(void *data, int j, int k)
  )
{
  if( L[L[0]] == 0 ) return; // singleton

  int p = 0; // previous node
  for(int i = R[0]; i; )
    {
    if( p == U[i] )
      {
      if( L[i] == 0 || L[L[i]] == 0 ) // skip to the right node
        {
        p = L[i];
        continue;
        }

      double S_LL_R = sim(data,L[L[i]],R[i]);
      double S_LR_R = sim(data,R[L[i]],R[i]);
      
      if( S_LL_R > S_LR_R )
        SWAP(int, L[L[i]], R[L[i]] );

      p = i;
      i = L[i];
      }
    else if( p == L[i] )
      {
      if( R[i] == 0 || L[R[i]] == 0 )
        {
        p = i;
        i = U[i];
        continue;
        }

      double S_L_RL = sim(data,L[i],L[R[i]]);
      double S_L_RR = sim(data,L[i],R[R[i]]);

      if( S_L_RL < S_L_RR )
        SWAP(int, L[R[i]], R[R[i]] );

      p = i;
      i = R[i];
      }
    else // if( p == R[i] )
      {
      p = i;
      i = U[i];
      }
    }
}

/*
 * convert nclust tree and scores to coordinates for line plots
 *
 */

void
dendrocoord (
  const int *root, // root of drawn subtree
  const int *U, const int *L, const int *R, const double *S,
  const int *leftmost, const int *nleaf,
  const int *nprune,    // branches with nleaf < *nprune are pruned
  const int *stemtype,  // 0: none, 1: relative, 2: absolute, 3: S[]
  const double *stemlength,
  double *x, double *y, int *npoints )  // max: 2*3*(N-1) + N
{
  int j = 0; // index to x and y coord   
  
  int p = U[*root]; // previous node
  for( int i = *root; i != U[*root]; )
    {
    if( p == U[i] ) // enter from parent
      {
      if( L[i] == 0                  // leaf node
          || nleaf[i] <= *nprune )   // pruned branch
        {                            // draw stem
        if( *stemtype == 1 )
          {
          x[j] = leftmost[i] + 0.5*nleaf[i]; y[j] = S[U[i]] + *stemlength; j++;
          }
        else if( *stemtype == 2 )
          {
          x[j] = leftmost[i] + 0.5*nleaf[i]; y[j] = *stemlength; j++;
          }
        else if( *stemtype == 3 )
          {
          x[j] = leftmost[i] + 0.5*nleaf[i]; y[j] = S[i]; j++;
          } 
        
        p = i;
        i = U[i];
        continue;
        } 

      x[j] = leftmost[i] + 0.5*nleaf[i]; y[j] = S[i]; j++;
      x[j] = leftmost[L[i]] + 0.5*nleaf[L[i]]; y[j] = S[i]; j++;
      p = i;
      i = L[i];
      }
    else if( p == L[i] )
      {
      x[j] = leftmost[L[i]] + 0.5*nleaf[L[i]]; y[j] = S[i]; j++;
      x[j] = leftmost[R[i]] + 0.5*nleaf[R[i]]; y[j] = S[i]; j++;
      p = i;
      i = R[i];
      }
    else // if( p == R[i] )
      {
      x[j] = leftmost[R[i]] + 0.5*nleaf[R[i]]; y[j] = S[i]; j++;
      x[j] = leftmost[i] + 0.5*nleaf[i]; y[j] = S[i]; j++;
      p = i;
      i = U[i];
      }
    }
  *npoints = j;
}

/*
 *  tree subset, modify in place
 */

void
subclustree (
  const int *z,        // subset indicator: if 0,  delete
  int *U,
  int *L,
  int *R,
  double *S,
  int *nleaf,
  int *order,
  int *leftmost,
  int *level,
  int *supleaf,
  int *Nnew_           // new number of leaves
  )
{
  //
  // First pass: flag deleted leaves and respective parents, relink to
  // siblings
  //
  int i, Nnew = 0;
  for(i = 1; L[i] == 0 ; i++ ) // loop over leaves (nodes with dummy children)
    {
    if( z[i-1] ) // length(z) == number of leaf nodes, index start at zero
      {
      Nnew++;
      continue;
      }
    int parent = U[i];
    int grandparent = U[parent];
    int sibling = (L[parent] == i ? R[parent] : L[parent]);
    int *C = (L[grandparent] == parent ? L : R );

    U[sibling] = grandparent;
    C[grandparent] = sibling;
    U[i] = U[parent] = -1;    // flag self and parent for deletion
    }

  //
  // Second pass: compact the nodes 
  //
  int N = i-1;
  int *new_id = (int*)calloc(2*N, sizeof(int));
  int j;
  for(j = 0, i = 0; i < 2*N; i++ )
    if( U[i] != -1 )
      new_id[i] = j++;
  
  // only for debug
  if( j != 2*Nnew ) 
    { mess("j != 2*Nnew\n"); goto END; }
  
  for(i = 0; i < 2*N; i++ )
    {
    if(U[i] == -1 ) continue;
    j = new_id[i];
    U[j] = new_id[U[i]];
    L[j] = new_id[L[i]];
    R[j] = new_id[R[i]]; 
    S[j] = S[i];
    }
  L[0] = R[0] = 2*Nnew - 1;
  
  branch_nleaf( U, L, R, nleaf );
  leafordering( U, L, R, order, leftmost );
  branch_level( &Nnew, U, level );

  for(int i = 1; i <= N; i++ )
    if( z[i-1] )
      supleaf[ new_id[i]-1] = i;

END:
  *Nnew_ = Nnew;
  free(new_id);
  return;
}

