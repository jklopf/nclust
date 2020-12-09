
#ifndef _NCLUSTREE_H_
#define _NCLUSTREE_H_

#ifndef SWAP
#define SWAP(type,x,y) do{ type t=(x);(x)=(y);(y)=t; } while(0)
#endif

/*
 * given a binary tree, returns the ordering of the leaves and
 * the left-most leaf of each cluster.
 */

extern void
leafordering( 
  const int *U,
  const int *L,
  const int *R,
  int *order,
  int *leftmost )
;

/*
 * count the number of leaves under each branch
 */
extern void
branch_nleaf (
  const int *U,
  const int *L,
  const int *R,
  int *nleaf        // the lenght is 2*N
  )
;

/*
 * level of the nodes, starting from 1 at the root
 *
 */
extern void
branch_level (
  const int *N_,
  const int *U,
  int *level
  )
;

extern void
branchflip_tightleft(
  int N,
  const int *U,
  int *L,
  int *R,
  const double *S
  )
;

extern void
branchflip_tightright(
  int N,
  const int *U,
  int *L,
  int *R,
  const double *S
  )
;

/*
 * nearest newphew branch flipping
 *
 */

extern void
branchflip_nnephew(
  const int *U,
  int *L,
  int *R,

  void *data,
  double (sim)(void *data, int j, int k)
  )
;

/*
 * dendrogram geometric coordinates from a tree
 */
void
dendrocoord (
  const int *root,
  const int *U, const int *L, const int *R, const double *S,
  const int *leftmost, const int *nleaf,
  const int *nprune,   // branches with < n leaves are pruned
  const int *stemtype,  // 0: none, 1: relative, 2: absolute, 3: S[]
  const double *stemlength,
  double *x, double *y, int *npoints )  // max: 2*3*(N-1) + N
;

/*
 * tree subset, modify in place
 *
 */
extern void
nclustree_subset (
  const int *z,        // subset indicator: if 0,  delete
  int *U,
  int *L,
  int *R,
  double *S,
  int *nleaf,
  int *order,
  int *leftmost,
  int *level,
  int *Nnew_           // new number of leaves
  )
;
#endif // _NCLUSTREE_H_
