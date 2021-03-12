ng_dendrogram <- function( clust, horizontal = FALSE, 
  wscale = NULL,  # scale of tree width (default 0.5..n+0.5)
  hscale = NULL,  # scale of tree height (range of dendrogram points)
  hexpand = 0.05, # expand by the given proportion
  hflip = FALSE,  # flip the direction of tree growth
  col = "black",
  ... )           # options for dendrocoord()
{
  u <- dendrocoord( clust, horizontal = horizontal, ... )
  nitems <- clust$N
  M <- cbind( c(1+hexpand,-hexpand),c(-hexpand,1+hexpand))

  if( horizontal )
    {
    if( is.null(wscale) ) wscale <- c(nitems+0.5, 0.5)
    if( is.null(hscale) ) hscale <- range(u$x)
    pushViewport( viewport(
      yscale = wscale,
      xscale = M %*% (if(hflip) rev(hscale) else hscale )) )
    }
  else
    {
    if( is.null(wscale) ) wscale <- c(0.5, nitems+0.5)
    if( is.null(hscale) ) hscale <- range(u$y)
    pushViewport( viewport(
      xscale = wscale,
      yscale = M %*% (if(hflip) rev(hscale) else hscale)) )
    }
  
  grid.clip()
  grid.lines(if(horizontal)u$x else u$x+0.5,
    if(horizontal) u$y+0.5 else u$y,default.units="native",
    gp=gpar(col=col))
  grid.clip(width=unit(100,"npc"),height=unit(100,"npc")) # ugly hack to undo clip
  return(hscale)
}
