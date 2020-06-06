#' Page faceting
#'
#' Wraper function for faceting page
#' @title page.facting
#' @param n.row 
#' @param n.col 
#' @param mar.lft 
#' @param mar.rht 
#' @param mar.top 
#' @param mar.btm 
#' @return faceting index
#' @author ZG Zhao
#' @export
page.faceting <- function(n.row, n.col, mar.lft=0, mar.rht=0, mar.top=0, mar.btm=0){
    if(missing(n.row) | missing(n.col)) stop('Please assign rows and columns.')

    xrow <- (100 - mar.top - mar.btm) %/% n.row
    lrow <- (100 - mar.top - mar.btm) %% n.row
    mar.btm <- mar.btm + lrow %/% 2
    mar.top <- mar.top + (lrow - lrow %/% 2)

    xcol <- (100 - mar.lft - mar.rht) %/% n.col
    lcol <- (100 - mar.lft - mar.rht) %% n.col
    mar.lft <- mar.lft + lcol %/% 2
    mar.rht <- mar.rht + (lcol - lcol %/% 2)

    mat <- matrix(rep(0, 10000), ncol=100)
    n <- 1
    
    for (i in 1:n.row){
        for (j in 1:n.col) {
            mat[(mar.top + 1 + (i - 1)*xrow):(mar.top + i * xrow),
                (mar.lft + 1 + (j - 1)*xcol):(mar.lft + j * xcol)] <- n
            n <- n + 1
        }
    }

    ndx.top <- 0
    if (mar.top > 0) {
        mat[1:mar.top, ] <- n
        ndx.top <- n
        n <- n + 1
    }

    ndx.btm <- 0
    if (mar.btm > 0) {
        mat[(100 - mar.btm + 1):100, ] <- n
        ndx.btm <- n
        n <- n + 1
    }
    
    ndx.left <- 0
    if (mar.lft > 0) {
        mat[(mar.top + 1):(100 - mar.btm), 1:mar.lft] <- n
        ndx.left <- n
        n <- n + 1
    }

    ndx.right <- 0
    if (mar.rht > 0) {
        mat[(mar.top + 1):(100 - mar.btm),(100 - mar.rht + 1):100] <- n
        ndx.right <- n
    }

    layout(mat)
    return(c(ndx.top, ndx.btm, ndx.left, ndx.right))
}

