\name{dendrocoord}
\alias{dendrocoord}
\title{Dendrogram coordinates of \verb{nclust} object}
\description{Convert the tree and clustering score in \verb{nclust}
  object into coordinates to be used in line plots}
\usage{
  dendrocoord( 
   clust,
   root = clust$R[1],
   height = "S",
   horizontal = FALSE,
   nprune = -400, stemtype = "none", stemlength = 0,
   f = identity )
}
\arguments{
  \item{clust}{\verb{nclust} object}
  \item{root}{Root of subtree to draw. The default is the whole
    tree.}
  \item{horizontal}{grow tree horizontally}
  \item{nprune}{Pruning of the tree. If positive,
    it is the maximum number of leafs on terminal clusters. If negative,
  it is the number of terminal clusters.}
  \item{stemtype}{The way the stemlength is interpreted: no stem is drawn
    for \verb{"none"}, \verb{"relative"} to the parent, at an
    \verb{"absolute"} (and constant) coordinate value, or at the position
    of the node's \verb{"score"}.  }
  \item{stemlength}{The length of the stem.}
  \item{f}{Function to transform \code{y} coordinates}
}
\details{
This function produces a list with members \verb{clust} and \verb{y}
that can be passed on to \code{lines(...)} or \code{plot(...,type="l")}
to draw the dendrogram.
}

