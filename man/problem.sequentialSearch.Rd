\name{problem.sequentialSearch}
\alias{problem.sequentialSearch}
\title{problem sequentialSearch}
\description{Compute the most likely peak model with at most the number of
peaks given by peaks.int. This function repeated calls
problem.PeakSegFPOP with different penalty values, until either
(1) it finds the peaks.int model, or (2) it concludes that there
is no peaks.int model, in which case it returns the next simplest
model (with fewer peaks than peaks.int).}
\usage{problem.sequentialSearch(problem.dir, peaks.int, verbose = 0)}
\arguments{
  \item{problem.dir}{problemID directory in which coverage.bedGraph has already been
computed. If there is a labels.bed file then the number of
incorrect labels will be computed in order to find the target
interval of minimal error penalty values.}
  \item{peaks.int}{int: target number of peaks.}
  \item{verbose}{Print messages?}
}

\value{Same result list from problem.PeakSegFPOP, with an additional
component "others" describing the other models that were computed
before finding the optimal model with peaks.int (or fewer)
peaks. Additional loss columns are as follows: under=number of
peaks in smaller model during binary search; over=number of peaks
in larger model during binary search; iteration=number of times
PeakSegFPOP has been run.}

\author{Toby Dylan Hocking}





