\name{rccenter}
\alias{ccenter,rcenter,rccenter}
\title{Row and/or column mean centering}
\description{Mean center the rows and/or the columns of a matrix, which
  can also be weighted.}
\usage{
  ccenter(x,w=NULL)
  rcenter(x,w=NULL)
  rccenter(x,w=NULL)
}
\arguments{
  \item{x}{Matrix to be centered}
  \item{w}{Weights in a matrix of the same dimension as \code{x}}
}
\details{
  \code{ccenter} centers the columns. \code{rcenter} centers the rows.
  \code{rccenter} centers both.
}
\value{
  Centered matrix.
}
