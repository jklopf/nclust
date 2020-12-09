wdownsample <- function(
  x,w=NULL,
  irow=NULL,      # row index (integer or characters)
  icol=NULL,      # col index (integer or characters)
  tilesize=c(1,1) # tile size
  )
{
  if( class(x)[1]=="tile" )
    {
    w <- data$w
    x <- data$x
    }
  else if(is.null(w)) 
    {
    return( downsample(x,irow,icol,tilesize))
    }

  mx <- nrow(x)
  nx <- ncol(x)
  if( is.null(irow) ) irow <- 1:mx
  if( is.null(icol) ) icol <- 1:nx

  if( !is.numeric(tilesize) )
    stop("tilesize has to be numeric")
  tilesize <- as.integer(tilesize)  
  if( any((tilesize < c(1,1)) | (tilesize > c(mx,nx))) )
    stop("tilesize out of range.")

  if( is.character(irow) )
    {
    if( is.null(rownames(x)) )
      stop("`x` has no rownames to match `irow`")
    irow <- match( irow, rownames(x) )
    }
  irow <- ifelse( is.na(irow) | irow > mx | irow < 1, 0, irow )

  # pad with null_id
  irow <- c(irow, rep(0, (tilesize[1] - length(irow)) %% tilesize[1]))

  if( is.character(icol) )
    {
    if( is.null(colnames(x)) )
      stop("`x` has no colnames to match `icol`")
    icol <- match( icol, colnames(x) )
    }
  icol <- ifelse( is.na(icol) | icol > nx | icol < 1, 0, icol )

  # pad with null_id
  icol <- c(icol, rep(0, (tilesize[2] - length(icol)) %% tilesize[2]))
  
  my <- length(irow) / tilesize[1]
  ny <- length(icol) / tilesize[2]

  r <- .C("wdownsample",
    as.integer(c(mx,nx,my,ny,tilesize)),
    as.integer(irow-1),
    as.integer(icol-1),
    as.double(x),
    as.double(w),
    y = double(my * ny),
    wy = double(my * ny),
    wrow = double(my),
    wcol = double(ny),
    NAOK=TRUE,
    PACKAGE="nclust")

  dim(r$y) <- c(my,ny)
  attr(r$y,"wrow") <- r$wrow
  attr(r$y,"wcol") <- r$wcol

  dim(r$wy) <- c(my,ny)
  attr(r$wy,"wrow") <- r$wrow
  attr(r$wy,"wcol") <- r$wcol

  list(y=r$y,wy=r$wy/outer(r$wrow,r$wcol))
}

wrow <- function(x) attr(x,"wrow")

wcol <- function(x) attr(x,"wcol")
