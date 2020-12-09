
# plot subset of labels, specified by regex and/or indices 
#
# - do not create new viewport, add labels to the outside of old ones!
#   This inherits the native coordinate of previous viewport.
# - 
#
ng_label1 <- function( 
  all,         # all labels (1:n)

  label,  # labels to be shown (a list of integer or character vectors)
               # if character, treat as regexp,  sandwiched between '^' and '$'
  col=NULL,     # vector of colors by the group names in the list
  horizontal=3, # 0 vertical, 1 horizonal, 2 tilted
  flip=FALSE,
  names = TRUE,
  wa = 1, wb = length(all),
  stem = unit(c(3,7,10,13,15),"pt"),
  base = unit(0,"pt")  # displacement from the default base
  )  
{
  n <- length(all)

  if(is.character(label) || is.integer(label))
    label <- list(label)  
  if(length(label) < 1) return()

  i <- lapply(label,function(g)
    if(is.character(g))
      unique(unlist( sapply(g, function(r) grep(paste0("^",r,"$"),all ))))
    else if(is.numeric(g))
      g <- g[ g >= 1 && g <= n ]
    else
      NULL
    )

  if(!is.null(col) )
    ci <- unlist(lapply(1:length(i),
      function(j) {
        cj <- col[ names(i)[j] ]
        if(is.null(names(i)) || is.na(cj)) cj <- "black"
        rep(cj,length(i[[j]]))
      }))
  else
    ci <- rep("black",sum(sapply(i,length)))

  i <- unlist(i)
  if(length(i) < 1) return()

  # clip to window
  ci <- ci[ i >= wa & i <=  wb ]
  i <- i[ i >= wa & i <= wb ]
  if(length(i) == 0 ) return()

  d <- duplicated(i)
  i <- i[!d]
  ci <- ci[!d]
  names(i) <- all[i]

  if(flip)
    { 
    base <- base + unit(0,"npc")   # bottom or left
    stempos <- base - rep( stem[-5], length(i) )
    labelpos <- base - rep( stem[5], length(i) )
    hjust <- 1
    }
  else
    {
    base <- base + unit(1,"npc")   # top or right
    stempos <- base + rep( stem[-5], length(i) )
    labelpos <- base + rep( stem[5], length(i) )
    hjust <- 0
    }

  grid.record(
      {
      ss <- 1.25 * c(1,1,1.414214,1.154701)[1+horizontal]
      j <- labelpack(
        i, 
        minx=ja, maxx=jb,
        space= ss * abs(as.numeric(convertUnit(
              unit(1,"strheight","X"), "native",
              typeFrom="dimension",
              axisFrom=ifelse(horizontal==0,"y","x")))))

      packstem <- unit( c(rbind(i,i,j,j)), "native" )

      if(horizontal == 0) 
        {
        yp <- packstem; xp <- stempos; 
        yt <- unit(j,"native"); xt <- labelpos
        }
      else
        {
        xp <- packstem; yp <- stempos; 
        xt <- unit(j,"native"); yt <- labelpos
        }

      grid.polyline(
        x=xp,
        y=yp,
        id.lengths=rep(4,length(i)),
        gp=gpar(col=col))

      grid.text( x=xt, y=yt, label=names(i),
        hjust=hjust, vjust=0.5, rot=c(0,90,45,60)[1+horizontal], gp=gpar(col=col))
      },

    list( i = i, col=ci, stempos = stempos, labelpos=labelpos, 
      ja=wa-0.5, jb=wb+0.5,horizontal=horizontal, 
      hjust= hjust)
    )

}

ng_label <- function(all,label,space=3,space.unit="cm",base=unit(0,"cm"),...)
{
  if(!is.list(label))
    {
    if( is.character(label) || is.numeric(label) )  # not a list
      label <- list(list(label))
    else
      stop("wrong label format (1)")
    }
  else if( all( sapply(label, function(u) is.character(u) || is.numeric(u))))
    {
    label <- list(label)
    }
  else if( ! all(sapply(label, 
        function(u) is.list(u) &&
          all(sapply(u,function(v) is.character(v) ||is.numeric(v))))))
    {
    stop("wrong label format (2)")
    }

  for(i in 1:length(label))
    ng_label1( all, label[[i]],base=base+unit((i-1)*space,space.unit),...)
}
