% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qPCR.R
\name{install.bioc}
\alias{install.bioc}
\title{Install Bioconductor packages}
\usage{
install.bioc(pkgs = character(), ...)
}
\arguments{
\item{pkgs}{character vector, names of packages to be installed}

\item{...}{params passing to \code{\link{BiocManager::install}}}
}
\description{
Install Bioconductor packages
}
\details{
"BiocManger" package will be installed if missing.
}
