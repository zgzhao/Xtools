#' A R warper for etags program
#' @title etags function
#' @param path
#' Path to source files
#' @param pattern
#' Source file subfix pattern
#' @param recursive
#' (TRUE, logical)
#' @param src
#' Vector of source file names
#' @param odir
#' Path to output TAGS file
#' @param append
#' (FALSE, logical)
#' @return NULL
#' @examples
#' require(xtools)
#' etags("/usr/local/share/emacs", pattern = "\\.el",
#'       odir="~/.emacs.d", append=FALSE)
#' @author ZGUANG
#' @export
etags <- function(path='.', pattern='\\.el$',
                  recursive=TRUE, 
                  src = list.files(path = path, pattern = pattern,
                       full.names = TRUE, recursive = recursive, ignore.case = TRUE),
                  odir = getwd(), append=FALSE){
    setwd(odir)
    if (append) {
        system(paste("etags -a", src[1]))
    } else {
        system(paste("etags", src[1]))
    }
           
    for (i in 2:length(src)) {
        system(paste("etags -a", src[i]))
    }
}

