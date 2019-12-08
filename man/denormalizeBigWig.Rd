\name{denormalizeBigWig}
\alias{denormalizeBigWig}
\title{denormalizeBigWig}
\description{Attempt to convert a bigWig file with non-negative non-integer
data to integer count data, by dividing each value by the smallest
non-zero value.}
\usage{denormalizeBigWig(input.bigWig, 
    output.bigWig, verbose = getOption("PeakSegPipeline.verbose", 
        1))}
\arguments{
  \item{input.bigWig}{Input bigWig file with non-integer data.}
  \item{output.bigWig}{Output bigWig file with integer data.}
  \item{verbose}{
}
}



\author{Toby Dylan Hocking}





