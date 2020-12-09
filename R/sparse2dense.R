sparse2dense <- function( p, transpose=FALSE )
{
  if(!is.data.frame(p) && ncol(p) < 2) 
    stop("input not a p list of data frame")
  
  p[[1]] <- as.character(p[[1]])
  p[[2]] <- as.character(p[[2]])
  
  if(transpose) p <- p[,c(2,1,3:ncol(p))]

  rnames <- unique(p[[1]])
  cnames <- unique(p[[2]])
  x <- array(0,dim=c(length(rnames),length(cnames)),
    dimnames=list(rnames,cnames))
  if(ncol(p)==2)
    for(i in 1:nrow(p))
      x[p[i,1],p[i,2]] <- 1
  else
    for(i in 1:nrow(p))
      x[p[i,1],p[i,2]] <- p[i,3]
  x
}

