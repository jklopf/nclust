pos2id <- function(clust,pos)
{
  if(is.character(pos))
    unlist(sapply(pos,function(s) grep(paste0("^",s,"$"),clust$labels)))
  else if(is.numeric(pos))
    clust$order[pos[which(pos > 0 & pos <= clust$N)]]
}


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
