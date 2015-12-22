#' Shiny app for qPCR data processing.
#'
#' Shiny app for qPCR data processing with BioConductor `HTqPCR` and `ddCt` packages.
#' @examples
#' library(Xtools)
#' qPCR.shiny.app()
#' @author ZG Zhao
#' @export
qPCR.shiny.app <- function(){
    a <- !require('HTqPCR', logical.return = T, quietly = T)
    b <- !require('ddCt', logical.return = T, quietly = T)
    if (a | b) {
        if (!require('BiocInstaller', logical.return = T, quietly = T)) {
            source("https://bioconductor.org/biocLite.R")
            biocLite()
        }
        if(a) biocLite('HTqPCR')
        if(b) biocLite('ddCt')
    }
                
    require(shiny)
    xpp <- system.file("examples", "qPCR", package="Xtools")
    runApp(xpp)
}

