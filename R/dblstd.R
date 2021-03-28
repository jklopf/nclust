dblstd <- function(y, iter_max = 10, tol=1e-5, verbose=0 )
{
  r <- .C("dblstd",
    dims=as.integer(c(nrow(y),ncol(y))),
    y=as.double(y),
    iter_max=as.integer(iter_max),
    tol=as.double(tol),
    verbose=as.integer(verbose)
    )
  dim(r$y) <- dim(y)
  dimnames(r$y) <- dimnames(y)
  attr(r$y,"n_iter") <- r$iter_max
  attr(r$y,"tol") <- r$tol
  r$y
}
