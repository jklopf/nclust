\name{dblstd}
\alias{dblstd}
\title{Standardize row and column}
\description{Scale the rows and columns simultaneously so that they
  all have the same root mean square}
\usage{
  dblstd(y, iter_max = 10, tol=1e-5, verbose=0)
}
\arguments{
  \item{y}{Data matrix}
  \item{iter_max}{Maximum number of iterations}
  \item{tol}{Tolerance (maximum deviantion from one of row RMS's)}
  \item{verbose}{If not zero, output iteration steps and score to stderr}
}
\details{
  Iteratively scale the rows and columns by dividing them by their respective
  root mean squares. Convergence is detected from the maximum deviation
  from one of the row RMS's.}
\value{
  Returns the scaled matrix.
}
