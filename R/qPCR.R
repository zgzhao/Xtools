#' Shiny app for qPCR data processing.
#'
#' Shiny app for qPCR data processing with BioConductor `HTqPCR` and `ddCt` packages.
#' @examples
#' library(Xtools)
#' qPCR.shiny.app()
#' @author ZG Zhao
#' @export
qPCR.shiny.app <- function(){
    a <- !library('HTqPCR', logical.return = T, quietly = T)
    b <- !library('ddCt', logical.return = T, quietly = T)
    if (a | b) {
        if (!library('BiocInstaller', logical.return = T, quietly = T)) {
            source("https://bioconductor.org/biocLite.R")
            biocLite()
        }
        if(a) biocLite('HTqPCR')
        if(b) biocLite('ddCt')
    }
                
    library(shiny)
    xpp <- system.file("examples", "qPCR", package="Xtools")
    runApp(xpp)
}

