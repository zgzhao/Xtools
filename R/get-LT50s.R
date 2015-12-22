.LT50sim <- function(TMP, EL, rt){
    EL[EL > 100] <- 100
    mT <- min(TMP) - 10
    
    ## Extending temperature low enough to ensure success logistic regression
    TMP <- c(TMP, mT - seq(0, by=5, length=5))
    EL <- c(EL, rep(100, 5))
    
    df <- na.omit(data.frame(TMP=TMP, EL=EL))
    
    mod0 <- nls(EL~ SSlogis(TMP, ASym, xmid, scal), data=df)
    R2 <- round(1 - var(residuals(mod0))/var(df$EL), 2)
    
    ## Predicted LT50 and 95% confidential intervals
    TMP.pred <- seq(rt[1], rt[2], by=0.01)
    ELs.pred <- predict(mod0, data.frame(TMP=TMP.pred))
    se.fit <- sqrt(apply(attr(ELs.pred, 'gradient'), 1, function(x) sum(vcov(mod0) * outer(x, x))))
    ELs.pred <- ELs.pred + outer(se.fit, qnorm(c(.5, .025,.975, 0.005, 0.995)))
    
    vals <- NULL
    for(i in 1:5) {
        els <- abs(ELs.pred[, i] - 50)
        xx <- TMP.pred[which(els==min(els))]
        vals <- c(vals, mean(xx))
    }
    
    pred <- cbind(TMP.pred, ELs.pred[, 1:3])
    colnames(pred) <- c('TMP', 'pred', 'low', 'high')
    result <- list(vals=vals, pred=pred)
    result
}


#' Calculate LT50s with electrolyte leakage data.
#'
#' User interface function to calculate LT50s with electrolyte leakage data. 
#' @title getLT50s() function
#' @param df Data frame of electrolyte leakage data (0~1). First column must be "TMP" for temperature.
#' @param range.t Temperature range for prediction.
#' @param all.info 
#' @return If all.info=FALSE (default), return LT50s and their CI95 only. Otherwise return a list consisting original data, LT50s and their CI95, and curve functions.
#' @author ZG Zhao
#' @export
getLT50s <- function(df, range.t=NULL, all.info=FALSE){
    TMP <- df$TMP
    if(is.null(range.t)) range.t <- range(TMP)
    vals <- NULL
    preds <- NULL
    for (i in 2:ncol(df)){
        EL <- df[, i] * 100
        result <- .LT50sim(TMP, EL, range.t)
        vals <- cbind(vals, result$vals)
        preds <- c(preds, list(result$pred))
    }
    
    names(preds) <- colnames(df)[-1]
    rownames(vals) <- c('LT50', 'CI95.L', 'CI95.U', 'CI99.L', 'CI99.U')
    colnames(vals) <- colnames(df)[-1]
    vals <- data.frame(vals)
    result <- list(data=df, vals=vals, preds=preds)
    class(result) <- 'LT50'
    
    if(all.info) return(result)
    else return(vals)
}



#' Calculate and plot LT50 curve.
#'
#' Calculate and plot LT50 curve.
#' @title Plot LT50 regression curves of given data.
#' @param dt The result of getLT50s() with 'all.info=TRUE'.
#' @param point.plot If TRUE (default), plot data points.
#' @param se.plot If TRUE (default), plot 95\% confidential intervals.
#' @param se.alpha Opacity of se shadows.
#' @param labels Replace sample names if set.
#' @param show.cols Sample columns to show in figure.
#' @param cex.point cex for data points.
#' @param cex.legend cex for legends.
#' @param xlab Same as 'plot' function.
#' @param ylab Same as 'plot' function.
#' @param xlim Same as 'plot' function.
#' @param ylim Same as 'plot' function.
#' @param lwd Same as 'plot' function.
#' @param lty Same as 'plot' function.
#' @param col Same as 'plot' function.
#' @param pch Same as 'plot' function.
#' @param ... Other parameters pass to 'plot' function.
#' @return NULL
#' @author ZG Zhao
#' @export
plot.LT50 <- function(dt, point.plot = TRUE, se.plot = TRUE, se.alpha=0.2, labels = NULL, show.cols=NA, 
                      cex.point=1, cex.legend=1,
                      xlab = expression(paste('Temperature (', degree, 'C)')), 
                      ylab = "Electrolyte Leakage (%)",
                      xlim = NULL, ylim = NULL, lwd=1, lty=NULL, col=NULL, pch=NULL, ...) {
    
    df <- dt$data
    preds <- dt$preds
    m <- length(preds)
    sels <- show.cols
    if(is.na(sels)) sels <- 1:m
    
    if(is.null(xlim)) xlim <- range(df$TMP) + c(-1, 1)
    if(is.null(labels) | length(labels)!=m) labels <- names(preds)
    cols <- col
    if(is.null(cols)) cols <- c('black', brewer.pal(m,"Dark2"))
    colx <- alpha(cols, se.alpha)
    pchs <- pch
    if(is.null(pchs)) pchs <- 1:m
    ltys <- lty
    if(is.null(ltys) | length(lty) != m) ltys <- rep(1, m)
       
    TMP <- seq(xlim[1], xlim[2], length=10)
    ELs <- seq(0, 110, length=10)
    
    plot(TMP, ELs, type='n', xaxs='i', xlim=xlim, xlab=xlab, ylab=ylab, ...)
    
    
    x <- df$TMP
    for(i in sels) {
        datax <- preds[[i]]
        tmps <- datax[, 1]
        els0 <- datax[, 2]
        points(tmps, els0, type='l', lwd = lwd, col=cols[i], lty=ltys[i])
        if(point.plot) {
            y <- unlist(df[, i + 1]) * 100
            y[y > 100] <- 100
            points(x, y, pch=pchs[i], col=cols[i], cex=cex.point)
        }
        if(se.plot) {
            els1 <- datax[, 3]
            els2 <- datax[, 4]
            polygon(c(tmps, rev(tmps)), c(els1, rev(els2)), border=colx[i], col=colx[i])
        }
    }
    
    legend('bottomleft', legend=labels[sels], pch=pchs[sels], lty=ltys[sels],
           col=cols[sels], inset=0.01, box.col=NA, cex=cex.legend)
}

#' isOverlap
#'
#' With two given intervals, return TRUE if they are overlap, otherwise return FALSE.
#' @title Is two intervals overlaped?
#' @param ir1 
#' @param ir2 
#' @return TRUE/FALSE.
#' @author ZG Zhao
#' @export
isOverlap <- function(ir1, ir2) {
    w1 <- max(ir1) - min(ir1)
    w2 <- max(ir2) - min(ir2)
    w0 <- max(c(ir1, ir2)) - min(c(ir1, ir2))
    w0 <= (w1 + w2)
}


#' Using LT50s tools with shiny.
#'
#' Using LT50s tools with shiny.
#' @examples
#' ## Begin example
#' library(Xtools)
#' LT50s.shiny.app()
#' ## End example
#' @author ZG Zhao
#' @export
LT50s.shiny.app <- function(){
    require(shiny)
    xpp <- system.file("examples", "LT50s", package="Xtools")
    runApp(xpp)
}
