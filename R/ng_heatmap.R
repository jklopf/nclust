require(grDevices)

ng_heatmap <- function(
  x, w=NULL,
  ro,   # row subset and ordering (indices or names)
  co,   # column subset and ordering (indices or names)
  rbin=1, cbin=10, # size of row and column bins for local averaging
  saturation=1,wsat=1,
  hcol=red.white.blue,
  rwb = nrow(x), cwa = 1
  )
{
  h <- wdownsample(x,w=w,irow=ro, icol=co, tilesize=c(rbin,cbin) )
  if( is.list(h) ) { y <- h$y; wy <- h$wy }
  else { y <- h; wy <- NULL }
  
  if(is.null(wy))
    {
    nc <- length(hcol)
    rh <- hcol[ nc %/% 2 + 1 + floor( (nc-1)/ 2 * tanh(saturation*y)) ]
    dim(rh) <- dim(y)
    }
  else
    {
    if(!is.matrix(hcol))
      hcol <- colormap2d("blue","gray","white","red")
    nh <- nrow(hcol)
    ns <- ncol(hcol)
    rh <- hcol[nh %/% 2 + 1 + floor( (nh-1)/2*tanh(saturation*y)) + nh*(
      1+floor((ns-1)*tanh(wsat*wy)))] 
    dim(rh) <- dim(y)
    }

  pushViewport(viewport(
      xscale=c(cwa - 0.5, cwa + length(co) + (cbin-length(co))%% cbin - 0.5 ),
      yscale=c(rwb + 0.5, rwb -( -0.5 + length(ro) + (rbin-length(ro)) %% rbin))
      ))

  grid.raster( rh,
    interpolate=FALSE,
    width=unit(1,"npc"),height=unit(1,"npc") )
}
