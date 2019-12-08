\name{problem.coverage}
\alias{problem.coverage}
\title{problem coverage}
\description{Ensure that coverage.bedGraph has been correctly computed for a
particular genomic segmentation problem.}
\usage{problem.coverage(problem.dir, 
    verbose = getOption("PeakSegPipeline.verbose", 
        1))}
\arguments{
  \item{problem.dir}{Path to a directory like sampleID/problems/problemID where
sampleID/coverage.bigWig contains counts of aligned reads in the
entire genome, and problemID is a chrom range string like
chr6_dbb_hap3:3491790-3736386 that indicates the genomic
coordinates of a particular segmentation problem. If
problemID/coverage.bedGraph does not exist, or its first/last
lines do not match the expected problemID, then we recreate it
from sampleID/coverage.bigWig.}
  \item{verbose}{print messages?}
}

\value{problem data.table. If necessary, the \code{\link{bigWigToBedGraph}} command
line program is used to create problemID/coverage.bedGraph and
then we (1) stop if there are any negative or non-integer data and
(2) add lines with zero counts for missing data.}

\author{Toby Dylan Hocking}





