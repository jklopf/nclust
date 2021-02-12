autocol <- function(k)
{
  if( k==1 ) return( c("black") )

  return( hsv( (0:(k-1))/k, s=.8,v=.8 ))
}


# color heatmaps of categorical labels

ng_tag <- function(
  z,  # tag matrix or data frame : row names should exist and match ref
  all,  # labels in the right order
  transpose=FALSE,  # column of z displayed as columns
  col.tag=NULL      # color scheme
  )
{
  if(transpose==FALSE)
    {
    m <- length(all)
    n <- ncol(z)
    }
  else
    {
    m <- ncol(z)
    n <- length(all)
    }

  if( !is.list(col.tag) ) # already a color matrix
    {
    h <- z[ match(all,rownames(z)), ]
    }
  else
    {
    h <- sapply( 1:ncol(z), function(j)
      {
        u <- z[,j]
        if(!is.factor(u)) u <- as.factor(u)
        k <- nlevels(u)

        if( is.null(col.tag) || is.null(col.tag[[colnames(z)[j]]]) )
          colj <- autocol(k)
        else
          colj <- col.tag[[colnames(z)[j]]]
        
        names(colj) <- as.character(levels(u))
        colj[ as.character(z[match(all,rownames(z)),j]) ]
      })
    }

  pushViewport(viewport())
  if(transpose==TRUE) h <- t(h)
  grid.raster( data.matrix(h), interpolate=FALSE, width=unit(1,"npc"),height=unit(1,"npc"))
  if( transpose )
    grid.text( y=(-0.5 + ncol(z):1)/(ncol(z)),
      x=unit(1,"npc")+unit(3,"pt"),label=names(z),rot=0,hjust=0,vjust=0.5)
  else
    grid.text( x=(-0.5 + 1:ncol(z))/(ncol(z)),
      y=unit(0,"npc")-unit(3,"pt"),label=names(z),rot=45,hjust=1,vjust=0.5)
}
