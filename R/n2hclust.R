# sort the order of the merge to match that produced by closest-pair
# algorithm. `stats::cutree` requires this
#
hclust_merge_order <- function(h)
{
  o <- order(h$height)
  h$height <- h$height[o]
  h$merge <- h$merge[o,]
  h$merge <- apply( h$merge,2,
    function(u) ifelse(u < 0,u,order(o)[ifelse(u<0,NA,u)]))
  h
}


# convert nclust output to 'hclust' format
#
n2hclust <- function( clust )
{
  N <- clust$N
  h <- list()
  class(h) <- "hclust"
  h$order <- clust$order
  o <- rev(order(clust$S[ (N+2):(2*N)]))
  h$height <- 1-clust$S[ o + N + 1 ]
  h$merge <- array( integer((N-1)*2),c(N-1,2))

  vid2hid <- c( seq( -1, -N, -1), rep( 0, N-1) )
  vid2hid[ N + o ] <- 1:(N-1)
  h$merge[,1] <- vid2hid[ clust$L[(N+2):(2*N)]][o]
  h$merge[,2] <- vid2hid[ clust$R[(N+2):(2*N)]][o]

  h$labels <- clust$labels

  hclust_merge_order(h)
}


# convert hclust tree to nclust
#
h2nclust <- function(hc)
{
  nc <- list()
  class(nc) <- "nclust"
  N <- dim(hc$merge)[1]+1
  nc$N <- N
  nc$U <- rep(0,2*N)
  nc$L <- rep(0,2*N)
  nc$R <- rep(0,2*N)
  nc$S <- rep(0,2*N)
  for(i in 1:(N-1))
    {
    nc$S[N+1+i] <- hc$height[i]
    Li <- hc$merge[i,1]
    Li <- ifelse(Li < 0, -Li, N+Li)
    nc$L[N+1+i] <- Li
    Ri <- hc$merge[i,2]
    Ri <- ifelse(Ri < 0, -Ri, N+Ri)
    nc$R[N+1+i] <- Ri 
    nc$U[Li+1] <- N+i
    nc$U[Ri+1] <- N+i
    }
  nc$L[1] <- nc$R[1] <- 2*N-1
  nc$labels <- hc$labels

  r <- .C("leafordering",
    as.integer(nc$U),
    as.integer(nc$L),
    as.integer(nc$R),
    order=integer(nc$N),
    leftmost=integer(nc$N*2),
    PACKAGE="nclust")
  nc$order <- r$order
  nc$leftmost <- r$leftmost

  r <- .C("branch_nleaf",
    as.integer(nc$U),
    as.integer(nc$L),
    as.integer(nc$R),
    nleaf=integer(2*N),
    PACKAGE="nclust")
  nc$nleaf <- r$nleaf
  
  r <- .C("branch_level",
    as.integer(N),
    as.integer(nc$U),
    level = integer(2*N),
    PACKAGE="nclust")
  nc$level <- r$level

  nc
}
