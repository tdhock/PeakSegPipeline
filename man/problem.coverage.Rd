\name{problem.coverage}
\alias{problem.coverage}
\title{problem coverage}
\description{Ensure that coverage.bedGraph has been correctly computed for a
particular genomic segmentation problem.}
\usage{problem.coverage(problem.dir)}
\arguments{
  \item{problem.dir}{Path to a directory like sampleID/problems/problemID where
sampleID/coverage.bedGraph contains counts of aligned reads in the
entire genome, and sampleID/problems/problemID/problem.bed
contains one line that indicates the genomic coordinates of a
particular segmentation problem.}
}

\value{Nothing. If necessary, the bigWigToBedGraph command line program
is used to create problemID/coverage.bedGraph and then we (1) stop
if there are any negative or non-integer data and (2) add lines
with zero counts for missing data.}

\author{Toby Dylan Hocking}





