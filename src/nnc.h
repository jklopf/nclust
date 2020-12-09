/*
 * nnc.h - Nearest-Neighbor Chain Algorithm for Hierarchical Clustering
 *
 */
#ifndef _NNC_H_
#define _NNC_H_

extern int    // 0: ok, -1: memory allocation failure
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
  int max_ncache,     // maximum cache length
  int verbose
  )
;

#endif // _NNC_H_
