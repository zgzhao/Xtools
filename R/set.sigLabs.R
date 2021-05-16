## get p-value of t.test for each pair of samples.
pairWisePval <- function(dls) {
    ## dls: list
    nn <- length(dls)
    snames <- names(dls)
    s1 <- s2 <- pp <- NULL
    for (i in 1:(nn - 1)) {
        dd1 <- dls[[i]]
        for (j in (i + 1):nn) {
            dd2 <- dls[[j]]
            s1 <- c(s1, snames[i])
            s2 <- c(s2, snames[j])
            if(length(dd1) != length(dd2))
                px <- t.test(dd1, dd2)$p.value
            else if(all(sort(dd1) == sort(dd2)))
                px <- 1
            else px <- t.test(dd1, dd2)$p.value
            pp <- c(pp, px)
        }
    }
    data.frame(S1 = s1, S2 = s2, P = pp)
}

## get p value from table for given samples.
getxp <- function(pTable, s1, s2) {
    ss <- unlist(apply(pTable[, 1:2], 1, FUN = function(x) length(setdiff(c(s1, s2), x)) == 0))
    rr <- pTable[[3]][ss]
    rr[1]
}

#' @export
getPmatrix <- function(pTable, sNames) {
    nn <- length(sNames)
    res <- NULL
    for(i in 1:nn) {
        pp <- NULL
        for(j in 1:nn)
            pp <- if (i == j) c(pp, 1) else c(pp, getxp(pTable, sNames[i], sNames[j]))
        res <- rbind(res, pp)
    }
    res
}

## set letter for samples
#' @export
setltts <- function(sNames, pTable){
    nsp <- getPmatrix(pTable, sNames) >= 0.05
    nsp <- as.list(as.data.frame(nsp))
    nn <- length(sNames)
    ltts <- rep("", nn)
    nx <- 1
    for(i in 1:nn) {
        if(ltts[i] != "") next
        aa <- nsp[[i]]
        ss <- sapply(nsp, FUN = function(bb) identical(aa, bb))
        ltts[ss & ltts == ""] <- letters[nx]
        nx <- nx + 1
    }
    ltts
}

## "signature" or fingerprint of a list
xsignature <- function(lst){
  nn <- length(lst)
  ft <- NULL
  for(i in 1:nn)
    for(j in 1:nn) {
      xss <- intersect(lst[[i]], lst[[j]])
      ft <- if(length(xss) < 1) c(ft, 0) else c(ft, 1)
    }
  paste(ft, collapse = "")
}

#' sigLabsPT: set significant labels for multiple comparisons with given p-value table
#'
#' Get/set significant labels for multiple comparisons. Samples with any identical letter are not significant different.
#' @title sigLabsPT
#' @param sNames sample names
#' @param pTable p-value of statistic test for each pair of samples.
#' @param scol1 column (name or index) of sample names 1
#' @param scol2 column (name or index) of sample names 2
#' @param pcol column (name or index) of of p-value.
#' @return string vector of significant letter
#' @author ZG Zhao
#' @export
#' @examples
#' library("Xtools")
#' pt <- read.csv(file.path(path.package("Xtools"), "data", "data.ptable.csv"))
#' str(pt)
#' dt <- subset(pt, subset=Gene=="NR"& Plant=="Arabidopsis")
#' smps <- sort(unique(c(dt$S1, dt$S2)))
#' cbind(sample=smps, siglabs=sigLabsPT(smps, dt, scol1=3, scol2=4, pcol=5))
sigLabsPT <- function(sNames, pTable, scol1 = 1, scol2 = 2, pcol = 3) {
    pTable <- pTable[, c(scol1, scol2, pcol)]
    if(! all(sNames %in% c(pTable[[1]], pTable[[2]])))
        stop("Sample names do not match the names in p-value table!")

    ## 1. set unique letter for each sample
    ltts <- setltts(sNames, pTable)
    lbs <- as.list(ltts)

    ## 2. add letter to sample with p<0.05 if not yet have an identical letter
    nsp <- getPmatrix(pTable, sNames) >= 0.05
    nsp <- as.list(as.data.frame(nsp))
    for(i in 1:length(ltts)) {
        lbx <- c(lbs[[i]], ltts[nsp[[i]]])
        lbs[[i]] <- sort(unique(lbx))
    }

    ## 3. remove redundant letters
    ## if "figureprint" does not change after removing the letter,
    ## then the letter is redundant.
    ft0 <- xsignature(lbs)
    lttx <- ltts
    for(aa in ltts) {
        if(! aa %in% lttx) next
        lbx <- lapply(lbs, FUN = function(x) setdiff(x, aa))
        if(any(sapply(lbx, FUN = function(x) length(x) < 1))) next
        if(xsignature(lbx) == ft0) {
            lbs <- lbx
            lttx <- setdiff(lttx, aa)
        }
    }
    ## 4. convert vector to string and adjust letter sequence
    lbs <- lapply(lbs, FUN = sort)
    lbs <- sapply(lbs, FUN = paste, collapse = "")
    if(length(ltts) > length(lttx)) {
        loo <- paste(lttx, collapse = "")
        lnn <- paste(letters[1:nchar(loo)], collapse = "")
        lbs <- chartr(loo, lnn, lbs)
    }

    lbs
}

#' Set significant labels of multi comparisons with given data.
#'
#' 1. get t-test p-value table of pairwise compairison; 2. set siglabs with function sigLabsPT.
#' @title sigLabsDT
#' @param dt data.frame. Support data of both long and short format.
#' @param col.data column indices or names for data. If only one column, the data.frame is treated as long format; otherwise the data.frame should be short format.
#' @param col.labels column indices or names for labels.
#' @param group.by column name or index for grouping. If set, siglabs will be set independently for each group.
#' @param stat.vals TRUE/FALSE(default). If FALSE (typically for short format data), return significant labels (vector) only; otherwise, return a data.frame.
#' @return string vector of significant letter
#' @author ZG Zhao
#' @export
#' @examples
#' library("Xtools")
#' dt <- read.csv(file.path(path.package("Xtools"), "data", "data.measure.csv"))
#' str(dt)
#' dx <- sigLabsDT(dt, col.data=4:7, group.by="Plant", stat.vals=TRUE)
#' dx
#' library("ggplot2")
#' px <- ggplot(dx, aes(x=Leaf, y=mean, fill=Treat)) + theme_bw(base_size = 16)
#' px + geom_bar(stat="identity", position=position_dodge(.9), width=0.8, color="gray20") +
#'   geom_errorbar(aes(ymin=mean - sd, ymax=mean + sd), position = position_dodge(.9), width=0.5) +
#'   geom_text(aes(y=mean + sd + offset, label=siglabs), position=position_dodge(.9), size=5) +
#'   scale_y_continuous(expand = c(0, 0)) + geom_blank(aes(y = (mean + sd) * 1.15)) +
#'   facet_wrap(~Plant, nrow=1)
sigLabsDT <- function(dt, col.data = 1:ncol(dt), col.labels=NULL,
                      group.by = NULL, stat.vals=FALSE) {
    cnames <- colnames(dt)
    if(is.numeric(col.data)) col.data <- cnames[col.data]
    if(is.null(col.labels)) {
        col.labels <- setdiff(cnames, col.data)
    } else if(is.numeric(col.labels)) {
        col.labels <- cnames[col.labels]
    }

    if(length(col.data) == 1) {
        ## long format data
        fml <- substitute(x ~ ., list(x=as.name(col.data)))
        dt <- aggregate(formula(fml), data=dt, FUN=c)
    } else {
        ## short format data
        dxx <- as.data.frame(t(dt[, col.data]))
        dxx <- as.list(dxx)
        names(dxx) <- NULL
        dt$CLEANDATA <- dxx
        col.data <- "CLEANDATA"
    }
    dtcheck <- lapply(dt[[col.data]], FUN=is.numeric)
    if(! all(unlist(dtcheck))) stop("Required numeric data columns!")

    if(is.null(group.by)) {
        dt$XTREATGROUP <- 1
        group.by <- "XTREATGROUP"
    } else if(is.numeric(group.by)) {
        group.by <- cnames[group.by[1]]
    } else group.by <- group.by[1]

    lbs <- rep("", nrow(dt))
    pts <- list()
    for(pp in unique(dt[[group.by]])) {
        ss <- dt[[group.by]] == pp
        dx <- dt[[col.data]][ss]
        xnames <- unlist(apply(dt[ss, col.labels], 1, paste0, collapse="-"))
        names(dx) <- xnames
        pt <- pairWisePval(dx)
        lbx <- sigLabsPT(xnames, pt)
        lbs[ss] <- lbx
        pts[[pp]] <- pt
    }
    if(stat.vals) {
        dt$mean <- unlist(lapply(dt[[col.data]], FUN=mean, na.rm=TRUE))
        dt$sd <- unlist(lapply(dt[[col.data]], FUN=sd, na.rm=TRUE))
        dt <- dt[, c(col.labels, "mean", "sd")]
        dt$siglabs <- lbs
        dt$offset <- 0
        gg <- dt[[group.by]]
        for(tt in unique(gg)){
            ss <- gg == tt
            dt$offset[ss] <- max(dt$mean[ss] + dt$sd[ss]) * 0.05
        }
        attr(dt, "pTables") <- pts
        dt
    } else lbs
}
