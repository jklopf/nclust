# flip the left-right ordering of a tree

flipclust <- function(C)
{
  C$order <- rev(C$order)
  if(!is.null(C$L))
    {
    L <- C$L; R <- C$R
    C$L <- R; C$R <- L
    }
  C
}
