# count the number of inversions immediately visible
#
n.inversion <- function(clust)
{
  i <- (clust$N+2):(2*clust$N-1)
  sum( clust$S[ clust$U[i]+1 ] > clust$S[i] )
}

# Check and fix inversions by lowering the score of the
# parent to match the child's.
# 
# Report the frequency of having to do this. the number is higher
# than that produced by `n.inversion` because a fix may introduce
# new violation.
#
check.inversion <- function(clust,verbose=1)
{
  N <- clust$N
  ninv <- 0
  maxd <- 0
  for(i in (N+2):(2*N-1))
    {
    d <- clust$S[ clust$U[i]+1 ] - clust$S[i]
    if( d > 0 )
      {
      if( d > maxd ) maxd <- d
      clust$S[ clust$U[i]+1 ] <- clust$S[i]
      ninv <- ninv + 1
      }
    }
  if( verbose && ninv )
    warning(paste0(ninv,"/",N-1,
        "(",sprintf("%.1f",100*ninv/(N-1)),
        "%) inversions found. max diff =",maxd,"."))
  clust$n.inversion <- ninv
  clust$maxdiff.inversion <- maxd
  invisible(clust)
}

list.inversion <- function(clust)
{
  N <- clust$N
  i <- which( (clust$S < clust$S[ clust$U+1 ])[(N+2):(2*N-1)] ) + N + 1
  Si <- clust$S[i]
  Sp <- clust$S[ clust$U[i]+1]
  cbind(
    node=i,
    S.node=Si,
    S.parent=Sp,
    abs.diff=abs(Si-Sp),
    "%.diff"=abs(Si-Sp)/abs(Si)
    )
}
