# summary of a dendogram

treesum <- function( C, label )
{
  o <- C$order
  S <- C$S[ C$U[ o + 1 ] + 1 ]
  b <- ifelse( C$L[ C$U[ o + 1 ] + 1] == o, "\\","/" )
  data.frame( item=label[o], score=S, branch=b, stringsAsFactors=FALSE )
}
