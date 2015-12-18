#* Function for error bar plotting
#' Plot error bars for scatter, line or barplots with arrows() function.
#' @title errorbar() function
#' @param x x values
#' @param y y values
#' @param sd.upr Standard errors for upper offset.
#' @param sd.lwr Standard errors for lower offset.
#' @param horiz Same as `horiz` in standard `barplot` function.
#' @param cex scale factor for blunt end
#' @param lwd Line width of error bar.
#' @param col Line color of error bar.
#' @return This fuction works as side effects.
#' @author ZG Zhao
#' @export
errorbar <- function(x, y, sd.upr, sd.lwr, horiz=FALSE, cex=1, lwd=1, col=1) 
{
    if(missing(sd.lwr) & missing(sd.upr)) return(NULL)
    if(all(sd.upr==0) & all(sd.lwr)) return(NULL)
    
    if(missing(sd.upr)) sd.upr <- sd.lwr
    if(missing(sd.lwr)) sd.lwr <- sd.upr
    
    if(!horiz){
        arrows(x, y, y1=y-sd.lwr, length=0.1*cex, angle=90, lwd=lwd, col=col)
        arrows(x, y, y1=y+sd.upr, length=0.1*cex, angle=90, lwd=lwd, col=col)
    } else{
        arrows(y, x, x1=y-sd.lwr, length=0.1*cex, angle=90, lwd=lwd, col=col)
        arrows(y, x, x1=y+sd.upr, length=0.1*cex, angle=90, lwd=lwd, col=col)
    }
}

