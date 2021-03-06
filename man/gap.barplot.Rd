% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gap.barplot.R
\name{gap.barplot}
\alias{gap.barplot}
\title{gap.barplot function}
\usage{
gap.barplot(
  df,
  y.cols = 1:ncol(df),
  sdu.cols = NULL,
  sdd.cols = NULL,
  btm = NULL,
  top = NULL,
  ylim = NULL,
  yext = 1.05,
  min.range = 10,
  max.fold = 5,
  ratio = 1,
  gap.width = 1,
  brk.side = 2,
  brk.type = "normal",
  brk.bg = "white",
  brk.srt = 135,
  brk.size = 1,
  brk.col = "black",
  brk.lwd = 1,
  error.cex = 1,
  error.lwd = 1,
  error.col = 1,
  sig.lab = NULL,
  sig.cex = 1,
  sig.xpos = 0.5,
  sig.ypos = -2,
  box.lwd = 1,
  box.col = 1,
  box.lty = 1,
  ...
)
}
\arguments{
\item{df}{long format of data.frame, grouped by columns}

\item{y.cols}{columns for y}

\item{sdu.cols}{}

\item{sdd.cols}{}

\item{btm}{bottom of gap}

\item{top}{top of gap}

\item{ylim}{}

\item{yext}{extension of y axis}

\item{min.range}{min ratio of max to min value, for automated calculation of gap}

\item{max.fold}{max ratio of max to bottom of gap, for automated calculation of gap.}

\item{ratio}{top to bottom spanning ratio}

\item{gap.width}{}

\item{brk.side}{if brk.size=1, set breaks on left only}

\item{brk.type}{normal or zigzag}

\item{brk.bg}{}

\item{brk.srt}{rotation angle of break ticks}

\item{brk.size}{}

\item{brk.col}{}

\item{brk.lwd}{}

\item{error.cex}{}

\item{error.lwd}{}

\item{error.col}{}

\item{sig.lab}{}

\item{sig.cex}{}

\item{sig.xpos}{}

\item{sig.ypos}{}

\item{box.lwd}{}

\item{box.col}{}

\item{box.lty}{}

\item{...}{to barplot function}
}
\value{
x value of bars
}
\description{
Bar plot with/without gaps
}
\details{
Set gap for barplot automatically or by btm and top pars.
}
\examples{
## Begin example
library(Xtools)
datax <- na.omit(airquality)[,1:4]
cols <- cm.colors(ncol(datax))
layout(matrix(1:6, ncol=2))
set.seed(0)
for (ndx in 1:6){
    dt <- datax[sample(rownames(datax), 10), ]
    par(mar=c(0.5,2,0.5,0.5))
    brkt <- sample(c('normal', 'zigzag'), 1)
    gap.barplot(dt, col=cols, brk.type=brkt, max.fold=5, ratio=2)
}
## End example
}
\author{
ZG Zhao
}
