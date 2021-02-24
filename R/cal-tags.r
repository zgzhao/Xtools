#* modified from shorten.to.string
short2str <-
    function(line, token)
{
    if (FALSE) {
        ans <- regexpr(strsplit("token", ",", fixed = TRUE)[[1L]][1L],
                       line, fixed = TRUE)
        if (ans == -1L) line
        else substr(line, 1L, ans + attr(ans, "match.length") - 1L)
    }
    else {
        ## can we just put essentially nothing? Seems to work
        substr(line, 1L, 1L)
    }
}

#* Modified from expr2token
stat.fun <-
    function(x,
             ok = c("<-", "=", "<<-", "assign", "setGeneric", "setGroupGeneric", "setMethod"),
             extended = TRUE)
{
    id <- ""
    value <-
        if ((length(x) > 1L) &&
            (length(token <- as.character(x[[2L]])) == 1L) &&
            (length(id <- as.character(x[[1L]])) == 1L) &&
            (id %in% ok)) token
        else
            character(0L)
    if (extended && identical(id, "setMethod"))
    {
        ## try to add the signature, comma separated
        sig <- tryCatch(eval(x[[3L]]), error = identity)
        if (!inherits(sig, "error") && is.character(sig))
            value <- paste(c(value, sig), collapse=",")
    }
    value
}

#* Modified from expr2token
stat.class <-
    function(x, ok = c("setClass", "setClassUnion"), extended = TRUE)
{
    id <- ""
    value <-
        if ((length(x) > 1L) &&
            (length(token <- as.character(x[[2L]])) == 1L) &&
            (length(id <- as.character(x[[1L]])) == 1L) &&
            (id %in% ok)) token
        else
            character(0L)
    value
}


#* Modified from rtags.file
cal.rtags.file <- function(src, type='fun')
{
    elist <- parse(src, srcfile = srcfile(src))
    if (length(elist) == 0) return(invisible(0))
    lines <- readLines(src)
    if(type=='fun') {tokens <- lapply(elist, stat.fun)} else {tokens <- lapply(elist, stat.class)}
    startlines <- sapply(attr(elist, "srcref"), "[", 1L)
    if (length(tokens) != length(startlines))
        stop("length mismatch: bug in code!", domain = NA)
    keep <- sapply(tokens, length) == 1L
    if (!any(keep)) return(invisible(0))
    return(length(tokens))
}

##* Modified from rtags
#' Calculate function or class number in source files
#' @title calTags function
#' @param path 
#' @param pattern 
#' @param recursive 
#' @param src 
#' @param cal.type String 'fun' or 'class'
#' @param keep.re 
#' @param verbose 
#' @return NULL
#' @author ZG Zhao
#' @export
calTags <-
    function(path = ".", pattern = "\\.[RrSs]$",
             recursive = FALSE,
             src = list.files(path = path,
             pattern = pattern,
             full.names = TRUE,
             recursive = recursive),
             cal.type='fun', 
             keep.re = NULL,
             verbose = getOption("verbose"))
{
    if (!missing(keep.re))
        src <- grep(keep.re, src, value = TRUE)
    n <- 0
    for (s in src)
    {
        if (verbose) message(gettextf("Processing file %s", s), domain = NA)
        tryCatch({n <- n + cal.rtags.file(s, cal.type)}, error = function(e) NULL)
    }
    invisible(n)
}

