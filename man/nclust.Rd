\name{nclust}
\alias{nclust}
\title{Hierarchical clustering}
\description{Hierarchical clustering using average-linked weighted dot-product}
\usage{
  nclust( y,
    w = NULL,
    wfeat=NULL,
    witem = NULL,
    method="avedot",
    branchflip = "center",
    standardize = TRUE,
    autofix.inversion = FALSE,
    cache_length = 32,
    verbose = 1 )
}
\arguments{
  \item{y}{Data matrix to be clustered. Rows are features and columns
    are items}
  \item{w}{Data matrix of weights, with the same dimensions as \option{y}.
  Default to all weights equal one.}
  \item{wfeat}{Marginal weights for features}
  \item{witem}{Marginal weights for items}
  \item{method}{\code{"avedot"} Average linkage of (optionally weighted)
    dot-product, or \code{"ward"} Ward linkage of Euclidean distance.}
  \item{branchflip}{Branch flipping methods: \code{"center"} (default),
    \code{"left"} or \code{"right"}.}
  \item{standardize}{Standardize the item vector, which turns covariance
  into (cosine) correlation similarity.}
  \item{autofix.inversion}{Fix inversion in tree heights.}
  \item{cache_length}{Maximum queue length for nearest neighbors. It does not
    affects the results, but only, possibly, computational efficiency.}
  \item{verbose}{0: quiet, +1: progress (default), +2 and +4: debugging}
}
\details{
Implement quadratic time and linear memory algorithm for hierarchical
clustering with weighted dot-product as similarity measure, with
average linkage.

If the data matrix is centered item-wise, the similarity measure
is the covariance. If it is also standardized, the similarity
measure is Pearson's correlation.

The binary tree output uses doubly-linked structure. Each node has a pointer
to the parent and both children. Although redundant, traversal algorithms
are simplified (e.g., no need for function recursion or iteration
with stack).
}
\value{
  A object of class \verb{nclust}, with members:
  \item{N}{Number of items}
  \item{L}{Index to the left node}
  \item{R}{Index to the right node}
  \item{U}{Index to the up (parent) node}
  \item{S}{Similarity score at each node}
  \item{order}{Ordering of items (indices to the original order)}
  \item{nleaf}{Number of leaves on each node}
  \item{leftmost}{The leftmost position of the leaf of a cluster}
  \item{level}{The level of a node (zero is the root)}
  \item{labels}{Leaf node labels (in input order)}
  \item{method}{Hierarchical clustering method}
  \item{branchflip}{Branch flipping method}
}
