##* modify from rtags.R of R base package
shorten.to.string <-
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

### * write.etags

## this function is responsible for formatting the output for a single
## file given the relevant information.  The format was inferred from
## the "Ctags" wikipedia entry and by studying etags output.

write.etags <-
    function(src,
             tokens, startlines, lines, nchars,
             ...,
             shorten.lines = c("token", "simple", "none"))
{
    ## extra 1 for newline
    shorten.lines <- match.arg(shorten.lines)
    offsets <- (cumsum(nchars + 1L) - (nchars + 1L))[startlines]
    lines <-
        switch(shorten.lines,
               none = lines,
               simple = sapply(strsplit(lines, "function", fixed = TRUE), "[", 1),
               token = mapply(shorten.to.string, lines, tokens))
    tag.lines <-
        paste(sprintf("%s\x7f%s\x01%d,%d",
                      lines, tokens, startlines,
                      as.integer(offsets)),
              collapse = "\n")
    ## simpler format: tag.lines <- paste(sprintf("%s\x7f%d,%d", lines, startlines, as.integer(offsets)), collapse = "\n")
    tagsize <- nchar(tag.lines, type = "bytes") + 1L
    cat("\x0c\n", src, ",", tagsize, "\n", tag.lines, "\n", sep = "", ...)
}

### * expr2token

## this computes the tag name from an expression.  Currently, this
## returns the second thing in the expression; so
##
##   foo <- function(x) ...      ==> `<-`,       foo, ...
##   setMethod("foo", "bar" ...  ==> setMethod,  foo, ...
##   setGeneric("foo", "bar" ... ==> setGeneric, foo, ...
##
## which covers the typical uses.  We match against a list to restrict
## types of expressions that are tagged.  To reject things like
##
##   x[i] <- 10
##
## the second component is required to have length 1.  One limitation
## is that things like
##
##   if (require(pkg)) foo <- ... else foo <- ...
##
## will not be handled.

expr2token <-
    function(x,
             ok = c("<-", "=", "<<-", "assign", "R6Class", 
                 "setGeneric", "setGroupGeneric", "setMethod",
                 "setClass", "setClassUnion"),
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

### * rtags.file

## Handles a single file

rtags.file <-
    function(src, ofile = "", append = FALSE,
             write.fun = write.etags) ## getOption("writeTags")
{
    
    ## FIXME: do we need to worry about encoding etc.?
    elist <- parse(src, srcfile = srcfile(src))
    if (length(elist) == 0) return(invisible())
    lines <- readLines(src)
    tokens <- lapply(elist, expr2token)
    startlines <- sapply(attr(elist, "srcref"), "[", 1L)
    if (length(tokens) != length(startlines))
        stop("length mismatch: bug in code!", domain = NA)
    keep <- lengths(tokens) == 1L
    if (!any(keep)) return(invisible())
    tokens <- unlist(tokens[keep])
    startlines <- startlines[keep]
    write.fun(src = src,
              tokens = tokens,
              startlines = startlines,
              lines = lines[startlines],
              nchars = nchar(lines, type = "bytes"),
              file = ofile, append = append)
}

### * rtags

#' Modified rtags function.
#'
#' Modified rtags function.
#' @title Modified rtags function.
#' @param path 
#' @param pattern 
#' @param recursive 
#' @param src 
#' @param keep.re 
#' @param ofile 
#' @param append 
#' @param verbose 
#' @return NULL
#' @author ZG Zhao
#' @export
xrtags <-
    function(path = ".", pattern = "\\.[RrSs]$",
             recursive = FALSE,
             src = list.files(path = path,
                              pattern = pattern,
                              full.names = TRUE,
                              recursive = recursive),
             keep.re = NULL,
             ofile = "", append = FALSE,
             verbose = getOption("verbose"))
{
    if (nzchar(ofile) && !append) {
        if (!file.create(ofile, showWarnings = FALSE))
            stop(gettextf("Could not create file %s, aborting", ofile),
                 domain = NA)
    }
    if (!missing(keep.re))
        src <- grep(keep.re, src, value = TRUE)
    for (s in src)
    {
        if (verbose) message(gettextf("Processing file %s", s), domain = NA)
        tryCatch(
                 rtags.file(s, ofile = ofile, append = TRUE),
                 error = function(e) NULL)
    }
    invisible()
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***

