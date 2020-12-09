tile <- function(...,transpose=NULL)
{
  ds <- list(...)
  n <- length(ds)
  names(ds) <- sapply(1:n,function(i) 
    ifelse(is.null(names(ds)[[i]]) || names(ds)[[i]]=="",
      paste0("#",i),names(ds)[[i]]))

  # check if dimnames is present for all
  for(i in 1:n)
    {
    if(is.null(dimnames(ds[[i]]))
      ||is.null(dimnames(ds[[i]])[[1]])
      ||is.null(dimnames(ds[[i]])[[2]])
      )
      stop(paste("Data set",names(ds)[[i]],"has incomplete dimnames"))
    }

  u1 <- c()
  u2 <- c()
  id1 <- list()
  id2 <- list()

  for(i in 1:n)
    {
    n1 <- dimnames(ds[[i]])[[1]]
    if( length(unique(n1)) != length(n1)) 
      stop(paste("dataset ",names(ds)[[i]]," dimnames[[1]] not unique."))
    n2 <- dimnames(ds[[i]])[[2]]
    if( length(unique(n2)) != length(n2)) 
      stop(paste("dataset ",names(ds)[[i]]," dimnames[[2]] not unique."))

    u1 <- union(u1,n1)
    u2 <- union(u2,n2)

    id1[[i]] <- match(n1,u1)
    id2[[i]] <- match(n2,u2)
    }

  names(id1) <- names(id2) <- names(ds)

  x <- array(0,dim=c(length(u1),length(u2)))
  w <- array(0,dim=c(length(u1),length(u2)))

  rownames(x) <- rownames(w) <- u1
  colnames(x) <- colnames(w) <- u2

  for(i in 1:n)
    {
    if(any( w[id1[[i]], id2[[i]]] != 0 ))
      stop(paste("Dataset ",names(ds)[[i]],": overlapping cells"))
    w[ id1[[i]], id2[[i]] ] <- ifelse(is.na(ds[[i]]),0,1)
    x[ id1[[i]], id2[[i]] ] <- ifelse(is.na(ds[[i]]),0,ds[[i]])
    }

  q <- list(x=x, w=w)
  if(all(w[1,1]==w)) q$w <- NULL
  class(q) <- "tile"
  invisible(q)
}
