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

\value{Either the data.table of coverage, or NULL if no coverage data exists.}

\author{Toby Dylan Hocking}





