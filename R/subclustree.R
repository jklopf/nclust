subclustree <- function(clust,include)
{
  if( length(include) == clust$N && is.logical(include) )
    {
    z <- include
    }
  else
    {
    if(is.character(include) && !is.null(clust$labels))
      {
      z <- rep(0,clust$N)
      z[ match(include,clust$labels) ] <- 1
      }
    else if(is.integer(include) && min(include) > 0 && max(include) <= clust$N )
      {
      z <- rep(0,clust$N)
      z[ include ] <- 1
      }
    else
      stop("can't understand the 'include'")
    }

  r <- .C("subclustree",
    z=as.integer(z),
    U=as.integer(clust$U),
    L=as.integer(clust$L),
    R=as.integer(clust$R),
    S=as.double(clust$S),
    nleaf=as.integer(clust$nleaf),
    order=as.integer(clust$order),
    leftmost=as.integer(clust$leftmost),
    level=as.integer(clust$level),
    supleaf=integer(length(clust$order)), # original leaf indices
    N=integer(1)
    ) 
  clust$N <- r$N
  NN <- 2*r$N
  clust$U <- r$U[1:NN]
  clust$L <- r$L[1:NN]
  clust$R <- r$R[1:NN]
  clust$S <- r$S[1:NN]
  clust$nleaf <- r$nleaf[1:NN]
  clust$order <- r$order[1:r$N]
  clust$leftmost <- r$leftmost[1:NN]
  clust$level <- r$level[1:NN]
  clust$supleaf <- r$supleaf[1:r$N]
  clust$labels <- clust$labels[ r$supleaf[1:r$N] ]
  invisible(clust)
}
