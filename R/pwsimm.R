wpack <- function(x,w)
{
  if(any(dim(x) != dim(w))) stop("dim(x) != dim(w)")
  wx <- array(c(w,x),dim=c(dim(x),2))
  aperm(wx,c(3,1,2))
}


pwsimm <- function(x,w=NULL)
{
  M <- nrow(x)
  N <- ncol(x)
  if(!is.null(w))
    {
    if(any(dim(w) != dim(x))) stop("dim(w) != dim(x)")
    W <- 2
    }
  else
    {
    if(length(dim(x)) == 2 ) W <- 1
    else if( length(dim(x)) == 3 && dim(x)[1]==2 )
      { W <- 3; M <- dim(x)[2]; N <- dim(x)[3] }
    else stop("wrong matrix format")
    }
  weighted <- ifelse(W > 1, 2,1)
  r <- .C("pwsimm",
    dims=as.integer(c(M,N,W)),
    as.double(x),
    as.double(w),
    S = double(N*(N+1)/2*weighted),
    PACKAGE="nclust"
    )  
  invisible(r$S)
}
