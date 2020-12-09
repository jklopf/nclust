

# convenience to navigate in a layout 
#
panel <- function(name, i,j)
{
  seekViewport(name)
  pushViewport(viewport(layout.pos.row=i,layout.pos.col=j))
}

