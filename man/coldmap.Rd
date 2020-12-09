\name{coldmap}
\alias{coldmap}
\title{Heatmap of large matrix}
\description{Cluster and show heatmap of a data matrix, with downsampling
  and dendrogram pruning}
\usage{
coldmap(
  x, w = NULL,
  rclust=NULL, cclust=NULL,
  clust=NULL,
  rwa=NULL,rwb=NULL,
  cwa=NULL,cwb=NULL,
  rlab=NULL,col.rlab=NULL,rlab.text=NULL,
  clab=NULL,col.clab=NULL,clab.text=NULL,
  rtag=NULL,col.rtag=NULL,
  ctag=NULL,col.ctag=NULL,
  tmarg = .5, bmarg=.5,
  lmarg = .5, rmarg=.5,
  rdend.space=5, cdend.space=4,
  rtag.space=0.5, ctag.space=0.5,
  rplot.space=0, cplot.space=1,
  rlab.space=3, clab.space=3,clab.tilt=FALSE,
	rdend.col = NULL,
  ltitle=NULL,
  ttitle=NULL,
  rstdz=TRUE,cstdz=TRUE,
  saturation = 2, wsat = 2, hcol=red.white.blue
  )
}
\arguments{
  \item{x}{Data matrix}
  \item{w}{Weight matrix}
  \item{rclust,cclust}{Row and column clustering objects, as produced
    by \code{nnc}}
  \item{clust}{A list with two members: \code{rclust} and \code{cclust}}
  \item{rwa,rwb,cwa,cwb}{Row and column windows boundaries measured in
    ordering indices. The automatic values are \code{rwa=1}, \code{rwb=nrow(x)},
  \code{cwa=1} and \code{cwb=ncol(x)}. Supplied values are clipped to
  be within these default ranges.}
  \item{rlab,clab}{Row and column labels, as specified for \code{ng_label}}
  \item{col.rlab,col.clab}{Color for \code{rlab} and \code{clab}}
  \item{rlab.text,clab.text}{Vector of labels. If not given, the row
    and column names of the matrix will be used. It is assumed to
  have the same ordering as the data matrix \code{x}}
  \item{rtag,ctag}{Annotations specified as data frames with color values.
    The columns correspond to variables, the rows  correspond to
    data rows (if \code{rtag}) or data columns (if \code{ctag}). The
  rows has to be named with corresponding data matrix rows or columns.}
  \item{tmarg,bmarg,lmarg,rmarg}{top, bottom, left, and right margins in
  centimeters}
  \item{rdend.space,cdend.space}{the space for row and column dendograms,
    in centimeters}
  \item{rdend.col,cdend.col}{Colors for the dendrograms}
  \item{ltitle,ttitle}{Left and top title, respectively}
  \item{rstdz,cstdz}{Use row and column standardization}
  \item{saturation}{Color saturation}
  \item{wsat}{Weight saturation}
  \item{hcol}{Heatmap color scheme}
}
\details{
  Given a data matrix, perform clustering of rows and columns. If either
  \code{rclust}, \code{cclust} and/or \code{clust} are specified, the
  supplied clustering objects are used instead. \code{clust} option
  overrides the other two. The function returns a list with
  \code{rclust} and \code{cclust} as members, and can be used
  in subsequent calls to this function.
}
\value{
  A list of two clustering objects (as produced by \code{nnc}), with
  two elements: \code{rclust} and \code{cclust}.
}
