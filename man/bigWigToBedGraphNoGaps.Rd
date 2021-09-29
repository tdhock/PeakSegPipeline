\name{bigWigToBedGraphNoGaps}
\alias{bigWigToBedGraphNoGaps}
\title{bigWigToBedGraphNoGaps}
\description{Run \code{\link{bigWigToBedGraph}} then fill in gaps with zeros.}
\usage{bigWigToBedGraphNoGaps(in.bigWig, 
    out.bedGraph, chrom = NULL, 
    start = NULL, end = NULL)}
\arguments{
  \item{in.bigWig}{character string path to input bigWig file.}
  \item{out.bedGraph}{character string path to output bedGraph file.}
  \item{chrom}{character string, chromosome name to filter data.}
  \item{start}{\code{start} position to filter data.}
  \item{end}{\code{end} position to filter data.}
}

\value{nothing.}

\author{Toby Dylan Hocking}





