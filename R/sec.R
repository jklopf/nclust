# convert to leaf node identifiers
# `pos` can be a character vector of labels or integer vector of dendrogram
#  coordinate
pos2id <- function(clust,pos)
{
  if(is.character(pos))
    unlist(sapply(pos,function(s) grep(paste0("^",s,"$"),clust$labels)))
  else if(is.numeric(pos))
    clust$order[pos[which(pos > 0 & pos <= clust$N)]]
}


# smallest enclosing cluster
# nodes are leaf node identifiers
# returning the branch identifier (0-based id)
#
sec <- function(
  clust,  # nclust object
  nodes   # node identifiers (2..N+1 for leafs, N+2..2N for branches)
  )
{
  commpath <- 2:(2*clust$N-1)
  for(i in nodes)
    {
      path <- c()
      j <- i
      while( j != 0 )
        {
        path <- c(path,j)
        j <- clust$U[j+1]
        }
      commpath <- intersect(commpath,path)
    }
  commpath[1]
}

# a set of (labeled) leaf indices under a smallest-encolsing cluster
# 
secitems <- function(
  clust,
  seeds,labels=TRUE
  )
{
  p <- pos2id( clust, seeds )
  if(length(p)==0) return(NULL)
  s <- 1+sec( clust, p )
  oa <- clust$leftmost[s]+1
  ob <- oa + clust$nleaf[s]-1;
  if(labels)
    {
    if( !is.null(clust$labels) )
      clust$labels[ clust$order[oa:ob]]
    else
      stop("cluster has no labels")
    }
  else
    clust$order[oa:ob]
}
