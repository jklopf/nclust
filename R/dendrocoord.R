dendrocoord <- function(
  clust,
  root = clust$R[1],
  height = "S",
  horizontal = FALSE,
  nprune = -400, stemtype="none", stemlength=0,
  f = identity )
{
  N <- clust$N
  if(root > clust$R[1] ) root <- clust$R[1]
  stemtype.i <-  pmatch(stemtype,c("none","relative","absolute","score"),nomatch=0)-1
  if( stemtype.i < 0 ) stop(paste0("unrecognized stemtype value: `",stemtype,"`"))
  if(height=="S") height <- clust$S
  else height <- log(1+clust$level)
  r <- .C("dendrocoord",
    root=as.integer(root),
    U=as.integer(clust$U),
    L=as.integer(clust$L),
    R=as.integer(clust$R),
    S=as.double(f(height)),
    leftmost=as.integer(clust$leftmost),
    nleaf=as.integer(clust$nleaf),
    nprune=as.integer(
      max(1,min(N,ifelse(nprune < 0,N/-nprune,nprune)))),
    stemtype=as.integer(stemtype.i),
    stemlength=as.double(stemlength),
    x=double(2*3*(N-1) + N),
    y=double(2*3*(N-1) + N),
    npoints=integer(1),
    PACKAGE="nclust")

  if( horizontal )
    xy <- list(y=r$x[1:r$npoints], x= r$y[1:r$npoints] )
  else
    xy <- list(x=r$x[1:r$npoints], y= r$y[1:r$npoints] )

  xy
}
