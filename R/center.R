# TODO: 
# - scale option (either row or column only)
# - factor for centering (rows and columns, but only on the margins)
#
ccenter <- function(x,w=NULL) 
{
  if(is.null(w)) 
    scale(x,scale=FALSE,center=TRUE)
  else
    {
    m <- apply(x*w,2,sum,na.rmt=T)/apply(w,2,sum,na.rm=T)
    scale(x,scale=F,center=ifelse(is.na(m),0,m))
    }
}

rcenter <- function(x,w=NULL)
{
  if(is.null(w))
    x - apply(x,1,mean,na.rm=T)
  else
    {
    m <- apply(x*w,1,sum,na.rm=T)/apply(w,1,sum,na.rm=T)
    x - ifelse(is.na(m),0,m)
    }
}

rccenter <- function(x,w=NULL) ccenter(rcenter(x,w=w),w=w)
