# obtain xcale from current viewport

ng_xscale <- function() as.numeric( convertX( unit(c(0,1),"npc"), "native") )
ng_yscale <- function() as.numeric( convertY( unit(c(0,1),"npc"), "native") )
