#' Read qPCR Ct data
#'
#' Ct values are exported from qPCR machine. Ct values should have labels of well index, gene name, sample tags and repeat number. Arrange of data:
#' - save data in csv or tab-delim format
#' - required columns: "Well", "Ct", "gene" and "rep.num"
#' - sample labels (required) can be put in a single "sample" column, or presented in multi columns named with "ds.xxxx" pattern.
#' - duplicated well label ("Well" column) are not allowed
#' - all samples should have tested same list of genes of identical repeats.
#' - avoid using space and other special character in column names.
#' @title read qPCR data table with Ct values. See "Details".
#' @param fname file name (path to Ct data file)
#' @param na.value Number. If Ct is NA or greater than this value, set to this value.
#' @param ... Parameters passed to read.table function.
#' @return list of data.ht (qPCRset) and data.ddct (ddCt) objects.
#' @author ZG Zhao
#' @export
read.qPCRtable <- function(fname, na.value=40, ...) {
    if(missing(fname)) stop("No input fname specified")
    dt.raw  <- read.table(fname, header = TRUE, ...)
    for(cc in colnames(dt.raw)) {
        dcc <- dt.raw[, cc]
        if (is.character(dcc)) dcc[is.na(dcc)] <- "NA"
        dt.raw[, cc] <- dcc
    }

    ## check column names
    cnames  <- colnames(dt.raw)
    if( sum(cnames == "Ct") < 1) stop("'Ct' column is required!")
    if( sum(cnames == "Well") < 1) stop("'Well' column is required!")
    if( sum(cnames == "gene") < 1) stop("'gene' column is required!")
    if( sum(cnames == "sample") < 1 && ! any(grepl("^ds\\.", cnames)))
        stop("Sample labels column(s) not found.")
    if( sum(cnames == "rp.num") < 1) stop("'rp.num' column for repeat number is required！")
    tbs <- table(dt.raw$rp.num)
    if(min(tbs) != max(tbs)) stop("Repeats should be identical accross all samples！")

    ## check samples and genes
    if (! any(grepl("^sample$", cnames))) {
        dt.raw$sample <- apply(dt.raw[, grep("^ds\\.", cnames)], 1, paste, collapse="@")
    }
    dt.raw$lb.sample <- paste(dt.raw$sample, dt.raw$rp.num, sep="@")
    tbs <- table(dt.raw$lb.sample)
    if(min(tbs) != max(tbs)) stop("All sample should tested identical list of genes!")
    tbs <- table(dt.raw$gene)
    if(min(tbs) != max(tbs)) stop("Each gene should be tested in all samples")

    # handle Ct values
    xct <- dt.raw$Ct
    xct <- gsub("(Undetermined)|(No Ct)", "60", xct, ignore.case = TRUE)
    xct <- as.numeric(xct)
    xct[xct > na.value] <- NA
    dt.raw$Ct <- xct

    genes <- unique(dt.raw$gene)
    dt.raw$gene <- factor(dt.raw$gene, levels=genes)
    dx <- dcast(dt.raw, gene~lb.sample,  value.var = "Ct")
    X <- as.matrix(dx[, -1])
    nspots <- nrow(X)
    nsamples <- ncol(X)

    samples <- sapply(colnames(X), function(x) {
        paste(rev(rev(strsplit(x, split="@", fixed = TRUE)[[1]])[ - 1]), collapse="@")
    })
    reps <- sapply(colnames(X), function(x) {
        rev(strsplit(x, split="@", fixed = TRUE)[[1]])[1]
    })

    ## qPCRset object =============================
    pdata <- data.frame(Sample=samples, rep=reps, row.names=1:nsamples)
    sample.info <- new("AnnotatedDataFrame", data = pdata)
    featPos <- paste0("A", 1:nspots)
    featType <- rep("Target", nspots)
    featType[grep("act|ref|ubq", genes, ignore.case = TRUE)] <- "Reference"
    featType <- factor(featType)

    df <- data.frame(featureNames=genes, featureType=featType, featurePos=featPos)
    metaData <- data.frame(labelDescription=c( "Name of the qPCR feature (gene)",
                                              "Type pf feature", "Position on assay"))
    featData <- AnnotatedDataFrame(data=df, varMetadata=metaData)

    colnames(X) <- NULL
    X.flags <- data.frame(matrix("Passed", ncol=nsamples, nrow=nspots))
    X.cat   <- data.frame(matrix("OK", ncol=nsamples, nrow=nspots),
                          stringsAsFactors=FALSE)
    htset <- new("qPCRset", exprs=X, phenoData=sample.info, featureData=featData,
        featureCategory=X.cat, flag=X.flags, CtHistory=data.frame())

    ## InputFrame object (ddCt) ====================
    coreData <- data.frame(matrix(nrow = nrow(dt.raw), ncol = 4))
    colnames(coreData) <- c("Sample", "Detector", "Ct", "Platename")

    ddset <- new("InputFrame")
    ddset@coreData <- coreData
    ddset@coreData$Sample <- dt.raw$sample
    ddset@coreData$Ct <- dt.raw$Ct
    ddset@coreData$Detector <- dt.raw$gene
    ddset@coreData$Platename <- rep(1, nrow(dt.raw))
    ddset@files <- fname

    list(data.ht=htset, data.ddct=ddset)
}

#' Read qPCR fluorescence data (Strategene qPCR machine only)
#'
#' Read qPCR fluo data exported by Strategene MxPro (Export Chart Data::To Text File::Format 1- Vertical grouped).
#' @title read fluorescence data
#' @param fname Character string: file name.
#' @param drop Character vector: Wells to drop, for example, drop=c("A1", "H") will drop A1 and all H wells.
#' @return data.frame.
#' @author ZG Zhao
#' @export
read.fluoData <- function(fname, drop=NULL) {
    df <- read.table(fname, skip = 2, sep="\t")[, 1:3]
    ndx <- grep("SYBR", df$V1)
    stt <- ndx + 1
    edd <- c(ndx[ - 1] - 1, nrow(df))
    pos <- gsub("^([a-z0-9]+).*", "\\1", df[ndx, 1], ignore.case = TRUE)
    for( i in 1:length(stt)) {
        xtt <- stt[i]
        xdd <- edd[i]
        df[xtt:xdd, 1] <- pos[i]
    }
    df <- df[ - ndx, ]
    colnames(df) <- c("Well", "Cycle", "Ct")
    for(cc in colnames(df)) {
        df[, cc] <- gsub("\\s+", "", df[, cc])
    }
    df$Cycle <- as.integer(df$Cycle)
    df$Ct <- as.numeric(df$Ct)
    well1 <- gsub("[^A-Z]", "", df$Well, ignore.case = TRUE)
    well2 <- gsub("[^0-9]", "", df$Well, ignore.case = TRUE)
    sel <- as.integer(well2) < 10
    well2[sel] <- paste0(0, well2[sel])
    df$Well <- paste0(well1, well2)

    df <- dcast(df, Cycle ~ Well, value.var = "Ct")
    if(! is.null(drop)) {
        cnames <- colnames(df)[ - 1]
        xcc <- toupper(drop)
        xcc <- paste(xcc, collapse="|", sep="")
        kep <- rep(TRUE, length(cnames))
        kep <- kep & !grepl(xcc, cnames)
        df <- df[, c(TRUE, kep)]
    }
    df
}

#' Calculate and return major results of HTqPCR.
#'
#' Return results includes d.norm, fold change and stats table.
#' @title Compare gene expression of sample pairs by limma method: FC and stats.
#' @param dt Data of qPCRset class, read with readDataHTqPCR function. Please refer to readDataHTqPCR function info.
#' @param ref.gene String. Names of reference gene(s).
#' @param comp.strings String of length 1. Compare stirng related to sample name(s).
#' @param comp.type Integer. 1: single control; 2: paired compares.
#' @param norm.method String. Default: del
#' @param verbose pass to \code{\link{HTqPCR::normalizeCtData}}
#' @return result: results of limmaCtData (list); contrast: contranst sample labels (named strings); d.norm: result of normalizeCtData.
#' @author ZG Zhao
#' @export
calHTqPCR <- function(dt, ref.gene, comp.strings=NULL, comp.type=1, norm.method="deltaCt", verbose=FALSE){
    results <- list(result=NULL, contrast=NULL, d.norm=NULL)
    if (is.null(dt)) return(results)
    res <- NULL
    vss <- NULL
    if ( ! ref.gene %in% fData(dt)$featureNames ) return(results)
    d.norm <- normalizeCtData(dt, deltaCt.genes=ref.gene, norm = norm.method, verbose=verbose)
    results$d.norm <- d.norm

    if (is.null(comp.strings)) return(results)
    samples    <- pData(dt)$Sample
    treatments <- gsub("[^a-z0-9]", ".", samples, ignore.case=TRUE)
    design     <- model.matrix(~0 + treatments)
    colnames(design) <- gsub('treatments', '', colnames(design))

    vxx <- gsub("(^\\s+|\\s+$)", "", comp.strings)
    if(vxx == "") return(results)
    if(comp.type==2) {
        if(!grepl(':', vxx)) return(results)
        sps <- strsplit(vxx, '[ ,;；，]+')[[1]]
        vxx <- sapply(sps, function(x){
            aa <- strsplit(x, ":", fixed = TRUE)[[1]]
            aa <- gsub("[^a-z0-9_]", ".", aa, ignore.case = TRUE)
        })
        if (all(vxx %in% treatments))
            vss <- apply(vxx, 2, paste, collapse=" - ")
        else return(results)
    } else {
        samples <- unique(samples)
        if(vxx %in% samples){
            sps <- samples[samples != vxx]
            trs <- gsub("[^a-z0-9_]", ".", sps, ignore.case = TRUE)
            vxx <- gsub("[^a-z0-9_]", ".", vxx, ignore.case = TRUE)
            vss <- paste(trs, vxx, sep=' - ')
        } else return(results)
    }
    names(vss) <- sps

    res <- NULL
    if(!is.null(vss)) {
        contrasts <- makeContrasts(contrasts=vss, levels=design)
        res <- limmaCtData(d.norm, design = design,
                           contrasts = contrasts, ndups = 1, spacing = 1)
    }
    results$result <- res
    results$contrast <- vss
    results
}

#' Install Bioconductor packages
#'
#' "BiocManger" package will be installed if missing.
#' @title Install Bioconductor packages
#' @param pkgs character vector, names of packages to be installed
#' @param ... params passing to \code{\link{BiocManager::install}}
#' @export
install.bioc <- function(pkgs=character(), ...){
    xpkgs <- .packages(all.available = TRUE)
    if(! "BiocManager" %in% xpkgs) install.packages("BiocManager")
    pkgs <- setdiff(pkgs, xpkgs)
    if(length(pkgs) > 0) BiocManager::install(pkgs, ask = FALSE, ...)
}

#' Update Bioconductor packages in your system
#'
#' Running this function will automatically upgrade the installed bioconductor packages in your system.
#' @title bioc.update
#' @return no return value
#' @author ZG Zhao
#' @export
bioc.upgrade <- function(){
    remove.packages(c("BiocManager", "BiocVersion"))
    install.packages("BiocManager")
    BiocManager::install(update=TRUE, ask=FALSE)
}

#' Shiny app for qPCR data processing.
#'
#' Shiny app for qPCR data processing with BioConductor `HTqPCR` and `ddCt` packages.
#' @examples
#' library(Xtools)
#' qPCR.shiny.app()
#' @author ZG Zhao
#' @export
qPCR.shiny.app <- function(){
    xpkgs <- .packages(all.available = TRUE)
    if(! all(c("HTqPCR", "ddCt") %in% xpkgs)) stop('Packages missing. Please install by running "install.bioc(c("HTqPCR", "ddCt"))"')
    library(shiny)
    xpp <- system.file("shiny", "qPCR", package="Xtools")
    runApp(xpp)
}
