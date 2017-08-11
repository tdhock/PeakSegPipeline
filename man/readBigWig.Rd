\name{readBigWig}
\alias{readBigWig}
\title{readBigWig}
\description{Read part of a bigWig file into R as a data.table (assumes
bigWigToBedGraph is present on your PATH).}
\usage{readBigWig(bigwig.file, chrom, start, end)}
\arguments{
  \item{bigwig.file}{path or URL of bigwig file.}
  \item{chrom}{chromosome to read.}
  \item{start}{position before reading.}
  \item{end}{plain text file where coverage is saved before reading into R.}
}

\value{data.table with columns chromStart chromEnd count.}

\author{Toby Dylan Hocking}





