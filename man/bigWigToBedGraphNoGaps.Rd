\name{bigWigToBedGraphNoGaps}
\alias{bigWigToBedGraphNoGaps}
\title{bigWigToBedGraphNoGaps}
\description{Run \code{\link{bigWigToBedGraph}} then fill in gaps with zeros.}
\usage{bigWigToBedGraphNoGaps(coverage.bigWig, 
    chrom, problemStart, 
    problemEnd, prob.cov.bedGraph)}
\arguments{
  \item{coverage.bigWig}{file/url of bigWig.}
  \item{chrom}{string, \code{chrom} of data to pull from bigWig.}
  \item{problemStart}{integer, start position of data to pull from bigWig.}
  \item{problemEnd}{integer, end position of data to pull from bigWig.}
  \item{prob.cov.bedGraph}{string, out file to save data.}
}

\value{nothing.}

\author{Toby Dylan Hocking}





