#ifndef _LWJ_H_
#define _LWJ_H_

extern int
lwj (
  int N,     // number of leaf nodes
  double *C, // [N*(N+1)/2] similarity matrix, diagonal unused

  int *L,     // [2*N] left
  int *R,     // [2*N] right
  int *U,     // [2*N] up (parent)
  double *S,  // [2*N] branch height

  int method, // 0: average
  int weighted, // 0: no, !=0: yes
  int max_ncache,
  int verbose
  );

#endif // _LWJ_H_
