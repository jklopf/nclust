\name{nclust}
\alias{nclust}
\title{Hierarchical clustering}
\description{Hierarchical clustering using average-linked weighted dot-product}
\usage{
  nclust( y, w = NULL, cache_length = 32, branchflip = 0,
    standardize = TRUE, verbose = 1 )
}
\arguments{
  \item{y}{Data matrix to be clustered. Rows are features and columns
    are items}
  \item{w}{Data matrix of weights, with the same dimensions as \option{y}.
  Default to all weights equal one.}
  \item{cache_length}{Maximum queue length for nearest neighbors. It does not
    affects the results, but only, possibly, computational efficiency.}
\item{branchflip}{Branch flipping methods. 0: nearest nephew (default), 1: tight-left, 2: tight-right}
  \item{standardize}{Standardize the item vector}
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
  \item{dim}{Array dimensions. Three dimensional if weighted, with
  the value-weight pair as the first dimensions.}
% \item{x}{Data array, including the weights if available}
  \item{L}{Index to the left node}
  \item{R}{Index to the right node}
  \item{U}{Index to the up (parent) node}
  \item{S}{Similarity score at each node}
  \item{order}{Ordering of items (indices to the original order)}
  \item{n}{Weights of each node}
  \item{leftmost}{The leftmost position of the leaf of a cluster}
  \item{level}{The level of a node}
  \item{options}{Options of the algorithm}
  \item{retstat}{Return status. O is OK}
}
