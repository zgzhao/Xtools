% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qPCR.R
\name{read.qPCRtable}
\alias{read.qPCRtable}
\title{read qPCR data table with Ct values. See "Details".}
\usage{
read.qPCRtable(fname, na.value = 40, ...)
}
\arguments{
\item{fname}{file name (path to Ct data file)}

\item{na.value}{Number. If Ct is NA or greater than this value, set to this value.}

\item{...}{Parameters passed to read.table function.}
}
\value{
list of data.ht (qPCRset) and data.ddct (ddCt) objects.
}
\description{
Read qPCR Ct data
}
\details{
Ct values are exported from qPCR machine. Ct values should have labels of well index, gene name, sample tags and repeat number. Arrange of data:
- save data in csv or tab-delim format
- required columns: "Well", "Ct", "gene" and "rep.num"
- sample labels (required) can be put in a single "sample" column, or presented in multi columns named with "ds.xxxx" pattern.
- duplicated well label ("Well" column) are not allowed
- all samples should have tested same list of genes of identical repeats.
- avoid using space and other special character in column names.
}
\author{
ZG Zhao
}
