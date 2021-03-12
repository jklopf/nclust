# heatmap of one matrix

ng_hm <- function(
  x, w = NULL,
  rclust=NULL, cclust=NULL,
  clust=NULL,
  method="avedot",
  rmethod=method,
  cmethod=method,
  
  rwa=NULL,rwb=NULL,
  cwa=NULL,cwb=NULL,

  rlab=NULL, clab=NULL,
  rlab.text=NULL, clab.text=NULL,
  rlab.space=3, clab.space=3,
  col.rlab=NULL, col.clab=NULL,
  clab.tilt=FALSE,

  rtag=NULL, ctag=NULL,
  col.rtag=NULL, col.ctag=NULL,
  rtag.space=0.5, ctag.space=0.5,

  tmarg = .5, bmarg=.5, lmarg = .5, rmarg=.5,
  
  rdend.space=5, cdend.space=4,
  rdend.col = NULL,cdend.col = NULL,
  
  rplot.space=0, cplot.space=1,

  ltitle=NULL, ttitle=NULL,
  rstdz=TRUE,cstdz=TRUE,
  saturation = 2, wsat = 2, hcol=red.white.blue
  )
{
  if(!is.matrix(x)) 
    stop("x is not a matrix")
  if( !is.numeric(x))
    stop("x is not numeric")

  # construct dendograms if not yet supplied
  if(!is.null(clust)) 
    { rclust <- clust$rclust; cclust <- clust$cclust }

  if(is.null(rclust))
    {
    if(is.null(w)) tw <- NULL
    else tw <- t(w)
    rclust <- nclust(t(x), w=tw, standardize=rstdz,method=rmethod )
    }
  
  if(is.null(cclust))
    cclust <- nclust(x, w=w, standardize=cstdz,method=cmethod )

  if(!is.null(rlab))
    {
    if(is.list(rlab) && is.list(rlab[[1]]))
      rmarg <- rmarg + rlab.space*length(rlab)
    else
      rmarg <- rmarg + rlab.space
    }

  if(!is.null(clab)) 
    {
    if(is.list(clab) && is.list(clab[[1]]))
      tmarg <- tmarg + clab.space*length(clab)
    else
      tmarg <- tmarg + clab.space
    }

  if(is.null(rtag))
    rtag.space <- 0

  if(is.null(ctag))
    ctag.space <- 0

  L <- grid.layout(
    ncol=6,
    widths=unit(c(lmarg, rdend.space,1,rtag.space,rplot.space,rmarg),
      c("cm","cm","null","cm","cm","cm")),
    nrow=6,
    heights=unit(c(cplot.space,ctag.space,tmarg,1,cdend.space,bmarg),
      c("cm","cm","cm","null","cm","cm"))
    )
  R_CTAG <- 2
  R_CLAB <- 3
  R_HM <- 4
  R_CDEND <- 5

  C_RDEND <- 2
  C_HM <- 3
  C_RTAG <- 4
  C_RLAB <- 5

  grid.newpage()
  pushViewport( viewport(layout=L, name="L"))

  # main coordinates: windows to cluster order
  m <- length(rclust$order)
  n <- length(cclust$order)
  if(is.null(rwa) || rwa < 1) rwa <- 1
  if(is.null(rwb) || rwb > m) rwb <- m
  if(is.null(cwa) || cwa < 1) cwa <- 1
  if(is.null(cwb) || cwb > n) cwb <- n

  # heatmap: with labels on the border
  panel("L",R_HM,C_HM)

  rbin = 1 + (rwb-rwa+1) %/% 500 
  cbin = 1 + (cwb-cwa+1) %/% 500
  
  if( !is.null(rclust$supleaf) ) ri <-rclust$supleaf
  else ri <- 1:rclust$N
  if( !is.null(cclust$supleaf) ) ci <- cclust$supleaf
  else ci <- 1:cclust$N

  ng_heatmap( x, w=w, 
    ri[ rclust$order[rwa:rwb]], ci[ cclust$order[cwa:cwb]],
    rbin=rbin,cbin=cbin,
    saturation = saturation, wsat = wsat, hcol=hcol, cwa = cwa, rwb = rwb )
  xs <- ng_xscale()
  ys <- ng_yscale()

  # row tags
  if( !is.null(rtag) )
    {
    panel("L", R_HM,C_RTAG)
    ng_tag( rtag, rownames(x)[ri[rclust$order[rwa:rwb]]], col.tag=col.rtag )

    #panel("L",R_HM-1,C_RTAG)
    #grid.text(x=( -0.5 + 1:ncol(rtag))/(ncol(rtag)),y=unit(3,"pt"),label=names(rtag),
    #  hjust=0,vjust=0.5,
    #  rot=90)
    }

  # row labels
  if( !is.null(rlab) )
    {
    panel("L", R_HM, C_RLAB)
    pushViewport(viewport(yscale=ys))
    if(!is.null(rlab.text)) all.rlab <- rlab.text[ri[rclust$order]]
    else if(!is.null(rownames(x))) all.rlab <- rownames(x)[ri[rclust$order]]
    else all.rlab <- c()
    ng_label( all.rlab, rlab, col=col.rlab, horizontal=0,
      wa=rwa,wb=rwb, base=unit(-1,"npc") )
    }
  
  # column tags
  if( !is.null(ctag) )
    {
    panel("L",R_CTAG,C_HM)
    ng_tag( ctag, transpose=TRUE,
      colnames(x)[ci[cclust$order[cwa:cwb]]], col.tag=col.ctag )
    #panel("L",R_CTAG,C_HM+1)
    #grid.text(x=unit(3,"pt"),y=(-0.5 + ncol(ctag):1)/(ncol(ctag)),
    #  label=names(ctag),hjust=0,vjust=0.5)
    }

  # column labels
  if( !is.null(clab) )
    {
    panel("L",R_CLAB,C_HM)
    pushViewport(viewport(xscale=xs))
    if(!is.null(clab.text)) all.clab <- clab.text[ci[cclust$order]]
    else if(!is.null(colnames(x))) all.clab <- colnames(x)[ci[cclust$order]]
    else all.clab <- c()
    ng_label( all.clab, space=clab.space, clab, col=col.clab,
      wa=cwa,wb=cwb, base=unit(-1,"npc"), horizontal=ifelse(clab.tilt,3,1) )
    }

  # row dendogram
  if( rdend.space > 0 && class(rclust)=="nclust" && !is.null(rclust$R) )
    {
    panel("L",R_HM,C_RDEND)
    pushViewport(plotViewport(c(0,4,0,1)))

		if( is.null(rdend.col) )
			rdend.col <- list( list(c(1,rclust$N), "black" ))
    else
      rdend.col <- c( list(list(c(1,rclust$N), "darkgray")), rdend.col )

    hrange <- range(rclust$S)
		lapply( rdend.col, function(u)
	      {
    		ng_dendrogram( rclust, horizontal=TRUE, wscale=ys,hscale=hrange,
     		stemtype="relative",stemlength=0.05, nprune=rbin,
				col = u[[2]], root = sec( rclust, pos2id(rclust,u[[1]])) 
				)
				})

    grid.xaxis(main=FALSE)
    grid.yaxis()
    }

  # column dendrogram
  if( cdend.space > 0 && class(cclust)=="nclust" && !is.null(cclust$R))
    {
    panel("L",R_CDEND,C_HM)
    pushViewport(plotViewport(c(2,0,1,0)))

    if( is.null(cdend.col) )
      cdend.col <- list( list(c(1,cclust$N), "black" ))
    else
      cdend.col <- c( list(list(c(1,cclust$N), "darkgray")), cdend.col )

    hrange <- range(cclust$S)
    lapply( cdend.col, function(u)
      {
      ng_dendrogram( cclust, horizontal=FALSE, wscale=xs,hscale=hrange,
        stemtype="relative",stemlength=0.05,nprune=cbin,
        col=u[[2]], root=sec(cclust,pos2id(cclust,u[[1]]))
      )
      })

    #ng_dendrogram( cclust, horizontal=FALSE, wscale=xs,
    # stemtype="relative",stemlength=0.05,nprune=cbin)
    
    grid.yaxis()
    grid.xaxis()
    }
  
  # margin titles

  if(!is.null(ltitle))
    {
    panel("L",R_HM,1)
    grid.text(ltitle,rot=90)
    }

  if(!is.null(ttitle))
    {
    panel("L",1,C_HM)
    grid.text(ttitle)
    }

  invisible( list(rclust=rclust,cclust=cclust) )
}

make_tag <- function( annot, varnames, cols=NULL )
{
  n <- length(varnames)
  if(is.null(cols)) cols <- rep("black",n)
  if(length(cols) < n)
    {
    temp <- cols
    cols <- rep("black",n)
    cols[1:length(temp)] <- temp
    }
  out <- do.call(cbind,lapply(1:n,function(i)
    data.frame(sapply(unique(annot[[varnames[i]]]),function(j)
        ifelse(annot[[varnames[i]]]==j,cols[i],NA)),stringsAsFactors=F)))
  rownames(out) <- rownames(annot)
  out
}
