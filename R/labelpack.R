labelpack <- function(x,space = NULL, spacefrac = 0.05,
  minx = min(x), maxx = max(x), maxiter = 50){
  n <- length(x)
  o <- order(x)

  if(missing(space))
    space <- spacefrac*(maxx-minx)
    
  r <- .C("R_labelpack",as.integer(n),as.double(x[o]),
    as.double(space),as.double(minx),as.double(maxx),
    as.integer(maxiter),y = double(n),
    PACKAGE="nclust")
  r$y[o] <- r$y
  r$y }

