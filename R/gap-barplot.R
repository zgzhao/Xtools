#' gap barplot
#'
#' gap barplot with various options such auto/manual gap interval, errorbar and significance labels.
#' @title gap.barplot function
#' @param df Long format data.frame, grouped by columns in barplot.
#' @param y.cols The index of y columns in df, such y.cols=1:3
#' @param sdu.cols If the df data frame contains standard error columns, sdu.cols assigns the column index.
#' @param sdd.cols If the upper and lower limits of SD are different, sdd.cols assigns the column index for lower limits.
#' @param btm Set the bottom break point value manually. Default is auto.
#' @param top Set the top break point value manually. Default is auto.
#' @param ylim y limits
#' @param yext If ylim is not set, yext sets the y extension.
#' @param min.range If the ratio of max y to min y is less than this value, the barplot will not have a gap.
#' @param max.fold The max ratio of max(all y) to max(lower part y).
#' @param ratio The height ratio of upper to lower part in barplot.
#' @param gap.width Relative width of the gap.
#' @param brk.type Break type, value should be `normal` or `zigzag`.
#' @param brk.bg Background color of the gap.
#' @param brk.srt Break line angle if brk.type is `normal`.
#' @param brk.size Break line length.
#' @param brk.col Break line color.
#' @param brk.lwd Break line width.
#' @param error.cex cex of error bar.
#' @param error.lwd lwd of error bar.
#' @param error.col color of error bar.
#' @param sig.lab Significance labels vector, the same order as y by columns.
#' @param sig.cex cex of sig labels.
#' @param sig.xpos Relative x adjustment of sig labels.
#' @param sig.ypos Relative y adjustment of sig labels.
#' @param box.lwd box line width
#' @param box.col box color.
#' @param box.lty box line type
#' @param ... Other parameters passed to barplot function.
#' @return Same as barplot return values.
#' @examples
#' ## Begin example
#' library(Xtools)
#' datax <- na.omit(airquality)[,1:4]
#' cols <- cm.colors(ncol(datax))
#' layout(matrix(1:6, ncol=2))
#' set.seed(0)
#' for (ndx in 1:6){
#'     dt <- datax[sample(rownames(datax), 10), ]
#'     par(mar=c(0.5,2,0.5,0.5))
#'     brkt <- sample(c('normal', 'zigzag'), 1)
#'     gap.barplot(dt, col=cols, brk.type=brkt, max.fold=5, ratio=2)
#' }
#' ## End example
#' @author ZG Zhao
#' @export
gap.barplot <- function(df, y.cols=1:ncol(df), sdu.cols=NULL, sdd.cols=NULL, btm=NULL, top=NULL, ylim=NULL, yext=1.05, 
                        min.range=10, max.fold=5, ratio=1, gap.width=1,  
                        brk.type='normal', brk.bg='white', brk.srt=135, brk.size=1, brk.col='black', brk.lwd=1,
                        error.cex=1, error.lwd=1, error.col=1,
                        sig.lab=NULL, sig.cex=1, sig.xpos=0.5, sig.ypos=-2, 
                        box.lwd=1, box.col=1, box.lty=1, ...){
    
    if (missing(df)) stop('No data provided.')

    ##=======================================================
    ## save options
    opts <- options()
    options(warn=-1)
    restore <- function()
    {
        options(opts)
        par(xpd=FALSE)
    }
    on.exit(restore())
    
    ##=======================================================
    ## set up data
    y <- t(df[, y.cols])
    colnames(y) <- NULL
    
    sdu <- y
    sdd <- y
    if (!is.null(sdu.cols)) sdu <- y + t(df[, sdu.cols])
    if (!is.null(sdd.cols)) sdd <- y - t(df[, sdd.cols])
    if (!is.null(sdu.cols) & is.null(sdd.cols)) sdd <- y - t(df[, sdu.cols])
    if (is.null(sdu.cols) & !is.null(sdd.cols)) sdu <- y - t(df[, sdd.cols])
    if(is.null(ylim)) {
        ylim <- c(0, max(sdu) * yext)
        if(any(y < 0)) ylim <- c(min(sdd) * yext, 0)
    }
    
    ##=======================================================
    ## break point1, apply barplot function if any y is negative.
    if(any(y < 0))
    {
        xx <- barplot(y, ylim=ylim, beside=TRUE, ...)
        if(!(is.null(sdu.cols) & is.null(sdd.cols))) errorbar(xx, y, sdu - y, y - sdd, cex=error.cex, lwd=error.lwd, col=error.col)   
        if(box.lwd > 0) box(lwd=box.lwd, col=box.col, lty=box.lty)
        if(!is.null(sig.lab))
            text(xx, sdu, labels=sig.lab, cex=sig.cex, adj=c(sig.xpos, sig.ypos), xpd=TRUE)
        return(invisible(xx))
    }
    
    ##=======================================================
    ## set up break point
    if (is.null(btm) | is.null(top)){
        autox <- .auto.breaks(sdu, min.range, max.fold)
        if (autox$flag){
            btm <- autox$btm
            top <- autox$top
        }
    }
    
    ##=======================================================
    ## break point2, apply barplot function if any y is negative.
    if (is.null(btm) | is.null(top))
    {
        xx <- barplot(y, ylim=ylim, beside=TRUE, ...)
        if(!(is.null(sdu.cols) & is.null(sdd.cols))) errorbar(xx, y, sdu - y, y - sdd, cex=error.cex, lwd=error.lwd, col=error.col)   
        if(box.lwd > 0) box(lwd=box.lwd, col=box.col, lty=box.lty)
        if(!is.null(sig.lab)) 
            text(xx, sdu, labels=sig.lab, cex=sig.cex, adj=c(sig.xpos, sig.ypos), xpd=TRUE)
        return(invisible(xx))
    }
        
    ## Set up virtual y limits
    halflen <- btm - ylim[1]
    xlen <- halflen * 0.1 * gap.width
    v_tps1 <- btm + xlen        # virtual top positions
    v_tps2 <- v_tps1 + halflen * ratio
    v_ylim <- c(ylim[1], v_tps2)
    r_tps1 <- top               # real top positions
    r_tps2 <- ylim[2]

    ## Rescale data
    lmx <- summary(lm(c(v_tps1, v_tps2)~c(r_tps1, r_tps2)))
    lmx <- lmx$coefficients

    sel1 <- y > top
    sel2 <- y >=btm & y <=top
    y[sel1] <- y[sel1] * lmx[2] + lmx[1]
    y[sel2] <- btm + xlen/2

    sel1 <- sdd > top
    sel2 <- sdd >=btm & sdd <=top
    sdd[sel1] <- sdd[sel1] * lmx[2] + lmx[1]
    sdd[sel2] <- btm + xlen/2

    sel1 <- sdu > top
    sel2 <- sdu >=btm & sdu <=top
    sdu[sel1] <- sdu[sel1] * lmx[2] + lmx[1]
    sdu[sel2] <- btm + xlen/2

    
    ## bar plot
    xx <- barplot(y, beside=TRUE, ylim=v_ylim, axes = FALSE, names.arg=NULL, ...)
    ## error bars
    errorbar(xx, y, sdu - y, y - sdd, cex=error.cex, lwd=error.lwd, col=error.col)

    ## Real ticks and labels    
    brks1 <- pretty(seq(0, btm, length=10), n=4)
    brks1 <- brks1[brks1 >= 0 & brks1 < btm]
    brks2 <- pretty(seq(top, r_tps2, length=10), n=4)
    brks2 <- brks2[brks2 > top & brks2 <= r_tps2]
    labx <- c(brks1, brks2)

    ## Virtual ticks
    brks <- c(brks1, brks2 * lmx[2] + lmx[1])
    axis(2, at=brks, labels=labx, ...)
    if(box.lwd > 0) box(lwd=box.lwd, lty=box.lty, col=box.col)
    
    ## break marks
    pos <- par("usr")
    xyratio <- (pos[2] - pos[1])/(pos[4] - pos[3])
    xlen <- (pos[2] - pos[1])/50 * brk.size
    px1 <- pos[1] - xlen
    px2 <- pos[1] + xlen
    px3 <- pos[2] - xlen
    px4 <- pos[2] + xlen
    py1 <- btm
    py2 <- v_tps1
    
    rect(px1, py1, px4, py2, col=brk.bg, border=brk.bg, xpd=TRUE)

    x1 <- c(px1, px1, px3, px3)
    x2 <- c(px2, px2, px4, px4)
    y1 <- c(py1, py2, py1, py2)
    y2 <- c(py1, py2, py1, py2)

    px <- .xy.adjust(x1, x2, y1, y2, xlen, xyratio, angle=brk.srt*pi/90)
    if (brk.type=='zigzag'){
        x1 <- c(x1, px1, px3)
        x2 <- c(x2, px2, px4)
        if (brk.srt > 90){
            y1 <- c(y1, py2, py2)
            y2 <- c(y2, py1, py1)
        } else {
            y1 <- c(y1, py1, py1)
            y2 <- c(y2, py2, py2)
        }
    }

    if (brk.type=='zigzag') {
        px$x1 <- c(pos[1], px2, px1, pos[2], px4, px3)
        px$x2 <- c(px2, px1, pos[1], px4, px3, pos[2])
        mm <- (v_tps1 - btm)/3
        px$y1 <- rep(c(v_tps1, v_tps1 - mm, v_tps1 - 2 * mm), 2)
        px$y2 <- rep(c(v_tps1 - mm, v_tps1 - 2 * mm, btm), 2)
    }
        
    par(xpd=TRUE)
    if(box.lwd > 0)
    {
        segments(px$x1, px$y1, px$x2, px$y2, lty=1, col=brk.col, lwd=brk.lwd)
    } else {
        nnn <- 1:length(px$x1)/2
        segments(px$x1[nnn], px$y1[nnn], px$x2[nnn], px$y2[nnn], lty=1, col=brk.col, lwd=brk.lwd)
    }
    
    ## significance labels
    if(!is.null(sig.lab)) text(xx, sdu, labels=sig.lab, cex=sig.cex, adj=c(sig.xpos, sig.ypos), xpd=TRUE)
    invisible(xx)
}

.xy.adjust <- function(x1, x2, y1, y2, xlen, xyratio, angle){
    xx1 <- x1 - xlen * cos(angle)
    yy1 <- y1 + xlen * sin(angle)/xyratio
    xx2 <- x2 + xlen * cos(angle)
    yy2 <- y2 - xlen * sin(angle)/xyratio
    return(list(x1=xx1, x2=xx2, y1=yy1, y2=yy2))
}

.auto.breaks <- function(dt, min.range, max.fold){
    ## 自动计算btm和top的思路
    ## 1. 判断是否需要设置断点的依据：最大值与最小值超过10倍
    ## 2. btm计算方法：将数据从小到大排序，从第2个数开始，如果它的值不超过全部数据最大值的5分之一，把它暂时设为btm，继续下一个直到满足前面条件。
    ## 3. top计算方法：初始值设为上方区域的第一个值的1/2,如果小于btm则分母和分子各加1,直到top>btm
    datax <- sort(unlist(dt))
    datax <- datax[datax != 0]
    n <- length(datax)
    flags <- FALSE
    btm <- top <- NULL
    if (length(datax) < 3) return(list(flag=flags, btm=btm, top=top))
    if (max(datax)/min(datax) < min.range) return(list(flag=flags, btm=btm, top=top))
    
    m <- max(datax)
    btm <- datax[2]
    i <- 3
    while(m/datax[i] > max.fold){
        btm <- datax[i]
        flags <- TRUE
        i <- i + 1
        if(i > n) break
    }

    if(i > n) i <- n
    
    if (flags) {
        btm <- btm + 0.05 * btm
        x <- 2
        top <- datax[i] * (x - 1)/x
        while (top < btm) {
            x <- x + 1
            top <- datax[i] * (x - 1)/x
            if (x > 100) {
                flags <- FALSE
                break
            }
        }
    }
            
    return(list(flag=flags, btm=btm, top=top))
}

