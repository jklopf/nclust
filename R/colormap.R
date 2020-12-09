colormap1d <- function(low,zero,high,n=10,...)
{
  c(
    colorRampPalette(c(low,zero),...)(n+1),
    (colorRampPalette(c(zero,high),...)(n+1))[-1]
  )
}

colormap2d <- function(low,missing,zero,high,ny=10,nw=20,...)
{
  z <- colorRampPalette(c(missing,zero),...)(nw)
  l <- colorRampPalette(c(missing,low),...)(nw)
  h <- colorRampPalette(c(missing,high),...)(nw)
  mapply(function(ll,zz,hh) colormap1d(ll,zz,hh,n=ny,...),l,z,h)
}

red.white.blue <- colormap1d("blue3","white","red3",n=10)

red.black.green <- colormap1d("green3","black","red3",n=10)



