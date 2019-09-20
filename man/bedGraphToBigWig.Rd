\name{bedGraphToBigWig}
\alias{bedGraphToBigWig}
\title{bedGraphToBigWig}
\description{Attempt to create \code{bigWig} file.}
\usage{bedGraphToBigWig(bedGraph, chromInfo, bigWig)}
\arguments{
  \item{bedGraph}{\code{bedGraph} file (input), four tab-separated columns: chrom,
chromStart, chromEnd, numeric data.}
  \item{chromInfo}{\code{chromInfo} file (input), two tab-separated columns: chrom, size in
bases.}
  \item{bigWig}{\code{bigWig} file (output).}
}

\value{TRUE if \code{bigWig} was created.}

\author{Toby Dylan Hocking}





