# TODO: 
# - tile object input and output
# - scale option (either row or column only)
# - factor for centering (rows and columns, but only on the margins), if
#   corresponding to a tile, better do it beforehand
# - c.f. ade4:::scalefacwt weights are marginal, and works only on columns
ccenter <- function(x,w=NULL) 
{
  if(is.null(w)) 
    scale(x,scale=FALSE,center=TRUE)
  else
    x - cbind(rep(1,nrow(x))) %*% (apply(x*w,2,sum,na.rm=T)/apply(w,2,sum,na.rm=T))
}

rcenter <- function(x,w=NULL)
{
  if(is.null(w))
    x - apply(x,1,mean,na.rm=T)
  else
    x - apply(x*w,1,sum,na.rm=T)/apply(w,1,sum,na.rm=T)
}

rccenter <- function(x,w=NULL) ccenter(rcenter(x,w=w),w=w)
