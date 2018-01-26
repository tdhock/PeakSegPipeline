\name{readCoverage}
\alias{readCoverage}
\title{readCoverage}
\description{Read sample coverage for one problem from either
sampleID/coverage.bigWig if it exists, or
sampleID/problems/problemID/coverage.bedGraph}
\usage{readCoverage(problem.dir, start, end)}
\arguments{
  \item{problem.dir}{project/samples/groupID/sampleID/problems/problemID}
  \item{start}{start of coverage to read.}
  \item{end}{end of coverage to read.}
}

\value{The data.table of coverage, which can be 1 row with 0 coverage if
no coverage data exists. (this is necessary for standard output of
the right number of samples for each joint problem, even those
with no coverage in some samples)}

\author{Toby Dylan Hocking}





