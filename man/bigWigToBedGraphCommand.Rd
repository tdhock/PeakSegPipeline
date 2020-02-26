\name{bigWigToBedGraphCommand}
\alias{bigWigToBedGraphCommand}
\title{bigWigToBedGraphCommand}
\description{Get command line to run \code{\link{bigWigToBedGraph}}.}
\usage{bigWigToBedGraphCommand(in.bigWig, 
    out.bedGraph, chrom = NULL, 
    start = NULL, end = NULL)}
\arguments{
  \item{in.bigWig}{character string path to input bigWig file.}
  \item{out.bedGraph}{character string path to output bedGraph file.}
  \item{chrom}{character string, chromosome name to filter data.}
  \item{start}{\code{start} position to filter data.}
  \item{end}{\code{end} position to filter data.}
}



\author{Toby Dylan Hocking}




\examples{
PeakSegPipeline::bigWigToBedGraphCommand(
  "in.bigWig", "out.bedGraph", "chr1", 0, 100)
}
