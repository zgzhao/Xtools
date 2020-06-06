#' Show pch values and symbols of plot function.
#'
#' Convenient function to show pch values and symbols of plot function.
#' @title pchShow
#' @return NULL
#' @author ZG Zhao
#' @export
pchShow <- function(){
    xx <- rep(1:10, 4)
    yy <- rep(4:1, each=10)
    pchs <- rep(1:20, 2)
    plot(xx, yy, pch=pchs, cex=2, col = (1:40 > 20) + 1, axes=FALSE, xlab='', ylab='')
    text(xx, yy+0.3, pchs, col="blue", xpd=TRUE, cex=0.8)
}
