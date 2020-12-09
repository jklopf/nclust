sim <- function(X,W=NULL)
{
  if(is.null(W)) weighted <- 0
  if(is.character(W) && W=="nna" ) weighted <- 1

  M <- nrow(X)
  N <- ncol(X)
  ws <- ifelse(weighted == 0, 1, 2)
  r <- .C("R_sim",
    as.integer(M),
    as.integer(N),
    as.double(X),
    C=double(N*(N+1)/2 * ws ),
    as.integer(weighted),
    double(0))
  attr(r$C,"N") <- N
  r$C
}
