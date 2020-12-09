\name{n2hclust}
\alias{n2hclust}
\title{Convert to \verb{hclust} object}
\description{Convert \verb{nclust} tree into \verb{hclust} tree
  format}
\usage{
  n2hclust( clust )
}
\arguments{
  \item{clust}{\option{nclust} object}
}
\details{
  The \verb{nclust} object produced by \verb{avdot} is converted
  (with loss of information) to an \verb{hclust} object amenable
  to various functions, such as dendrogram drawing and cutting.}

