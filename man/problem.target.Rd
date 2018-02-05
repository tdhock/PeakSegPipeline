\name{problem.target}
\alias{problem.target}
\title{problem target}
\description{Compute target interval for a segmentation problem. This function
repeated calls problem.PeakSegFPOP with different penalty values,
until it finds an interval of penalty values with minimal label
error. The calls to PeakSegFPOP are parallelized using mclapply if
you set options(mc.cores).}
\usage{problem.target(problem.dir, minutes.limit = getOption("PeakSegPipeline.problem.target.minutes", 
    Inf))}
\arguments{
  \item{problem.dir}{problemID directory in which coverage.bedGraph has already been
computed. If there is a labels.bed file then the number of
incorrect labels will be computed in order to find the target
interval of minimal error penalty values.}
  \item{minutes.limit}{Time limit; the search will stop at a sub-optimal target interval
if this many minutes has elapsed. Useful for testing environments
with build time limits (travis).}
}

\value{List of info related to target interval computation: target is the
interval of log(penalty) values that achieve minimum incorrect
labels (numeric vector of length 2), target.iterations is a
data.table with target intervals as a function of iteration,
models is a data.table with one row per model for which the label
error was computed.}

\author{Toby Dylan Hocking}





