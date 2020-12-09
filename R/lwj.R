lwj <- function(X,W=NULL,
  flip="nneph",method="ave",verbosity=1)
{
  N <- ncol(X)
  M <- nrow(X)
  
  i <- pmatch(flip,c("merge","tightleft","tightright","nnephew"))
  if(is.na(i)) flip <- 1
  else flip <- i-1
  
  i <- pmatch(method,c("ave","ward","single","complete"))
  if(is.na(i)) method <- 0
  else method <- i-1

  if( !is.null(W) )
    {
    if( W == "nna" ) { weighted <- 1; W <- double(0) }
    else if( is.matrix(W) && all(dim(W) == dim(X)) )
      weighted <- 2
    else stop("improper weight")
    }
  else { weighted <- 0; W <- double(0) }


  r <- .C("R_lwj",
    M = as.integer(M),
    N = as.integer(N),
    weighted = as.integer(weighted),
    X = as.double(X),
    W = as.double(W),
    L = integer(2*N),
    R = integer(2*N),
    U = integer(2*N),
    S = double(2*N),
    order = integer(N),
    nleaf = integer(2*N),
    leftmost = integer(2*N),
    level = integer(2*N),
    method = as.integer(method),
    flip = as.integer(flip),
    verbosity = as.integer(verbosity),
    NAOK=TRUE,
    PACKAGE="nclust"
    )
  r$X <- NULL
  r$labels <- colnames(X)
  class(r) <- "nclust"
  r
}
