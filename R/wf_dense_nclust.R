wf_dense_nclust <- function(
  data,
  w = NULL,
  wfeat = NULL, # feature weights
  witem = NULL, # item weights
  method = "average",
  branchflip = "center",
  standardize = TRUE,
  autofix.inversion = FALSE,
  cache_length = 32,
  verbose = 1
  )
{
  x <- NULL
  if( class(data)[1] == "tile" )
    {
    x <- data$x
    if(!is.null(data$w)) w <- data$w
    }
  else if(is.matrix(data)) 
    x <- data

  if( length(dim(x)) != 2 ) stop("length(dim(x)) != 2")
  if( !is.null(w) && any(dim(w) != dim(x)) ) stop("dim(w) != dim(x)")
  M <- nrow(x)
  N <- ncol(x)
  if( cache_length < 1 ) cache_length <- 1

  xx <- cbind( rep(0,M), x, array( 0, dim=c(M,N-1) ))
  if( is.null(w) )
    ww <- NULL
  else
    ww <- cbind( rep(0,M), w, array( 0, dim=c(M,N-1) ))
  
  if(is.null(wfeat))
    tt <- rep(1,M)
  else
    {
    if(length(wfeat) != M) stop("length(wfeat) != M")
    tt <- ifelse(wfeat < 0, 0, wfeat)
    }

  if(is.null(witem))
    uu <- c(0,rep(1,N),rep(0,N-1))
  else
    {
    if(length(witem) != N ) stop("length(witem) != N")
    uu <- c(0, ifelse(witem < 0, 0, witem), rep(0,N-1))
    }

  if(is.na(method <- pmatch(method,c("average","ward"))))
    stop("invalid `method`")

  if( method == 2 && !is.null(w) )
    stop("weighted data not implemented for Ward's method")

  if(is.na(branchflip <- pmatch(branchflip,c("center","left","right"))))
    stop("invalid `branchflip`")

  r <- .C("wf_dense_nclust",
    dims = as.integer( c(ifelse(is.null(w),1,2),M,N) ),
    options = as.integer( 
      c(cache_length, branchflip, standardize, verbose, method) ),

    ## input
    xx = as.double(xx),
    ww = as.double(ww), # ?
    tt = as.double(tt),
    uu = as.double(uu),

    ## tree output
    L = integer(2*N), R = integer(2*N), U = integer(2*N),
    S = double(2*N),
    order = integer(N),
    nleaf = integer(2*N),
    leftmost = integer(2*N),
    level = integer(2*N),

    retstat = integer(1), # return status
    NAOK=FALSE,
    PACKAGE="nclust"
    ) 
  if(r$retstat != 0 ) stop("nnc failed")
  r$retstat <- NULL

  r$dims <- r$dims[3]
  names(r)[1] <- "N"
  r$options <- NULL
  r$xx <- NULL
  r$ww <- NULL
  r$tt <- NULL
  r$uu <- NULL
  if( !is.null(colnames(x)) ) r$labels <- colnames(x)
  r$method <- c("average","ward")[method]
  r$branchflip <- c("center","left","right")[branchflip]

  rr <- check.inversion(r,verbose=verbose)
  if(autofix.inversion && rr$n.inversion > 0) r <- rr
  class(r) <- "nclust"
  invisible(r)
}
