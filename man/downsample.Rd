\name{downsample}
\alias{downsample}
\alias{wrow}
\alias{wcol}
\title{Means of a tiled subset of a matrix}
\description{Subset a matrix by arbitrarily ordered indices and
  calculate the means of its regular tiles}
\usage{
  downsample(x, irow=NULL, icol=NULL, tilesize=c(1,1) )
  wrow(x)
  wcol(x)
}
\arguments{
  \item{x}{input matrix}
  \item{irow}{indices to rows (integers or names)}
  \item{icol}{indices to columns (integer or names)}
  \item{tilesize}{tile row and column size}
}
\details{
  A subset of matrix \var{x} is created based on the indices,
  which may be contains arbitrary re-ordering and replicate
  references, as well as missing references (such as out-of-range
  integers or unmatched names). A missing reference produces NA's
  for the entire row or column.

  The subset is then ``tiled'', or divided into regular blocks
  with the given tile size. Each tile is averaged, taking into
  account the missing reference and possible missing values
  in \var{x}. The last tile in each row or column is padded by NA.

  The number of non-missing row and column references in each tile 
  is returned as attribute \var{wrow} and \var{wcol}, that can be
  conveniently accessed by the function \code{wrow(x)[i]} 
  and \code{wcol(x)[j]}.
  If there is no missing values in \var{x}, the weight of the cell \verb{y[i,j]}
  is \verb{wrow(y)[i] * wcol(y)[j]}. The weight in the presence of
  NA in the cells of \verb{x} is not kept. It is the number of
  non-missing values in each tile before taking the mean.
}
\value{
  A matrix of tile means, with attribute \var{wrow} and \var{wcol}
  corresponding to vectors of the number of non-missing references.
  }
\examples{
  x <- cbind( c(1,2), c(3,4) )
  downsample(x) # pass through
  downsample(x,tilesize=c(2,2)) # return 1x1 matrix of the grand mean
  downsample(x,irow=c(1,2,1,2),icol=c(2,2,1,1))  # duplicate, re-order
  downsample(x,irow=c(1,2,1),icol=c(2,1),tilesize=c(2,2))
}
