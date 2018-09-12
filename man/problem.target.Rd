\name{problem.target}
\alias{problem.target}
\title{problem target}
\description{Compute target interval for a segmentation problem. This function
repeated calls problem.PeakSegFPOP with different penalty values,
until it finds an interval of penalty values with minimal label
error. The calls to PeakSegFPOP are parallelized using mclapply if
you set options(mc.cores).}
\usage{problem.target(problem.dir, minutes.limit = NULL, verbose = 0)}
\arguments{
  \item{problem.dir}{problemID directory in which coverage.bedGraph has already been
computed. If there is a labels.bed file then the number of
incorrect labels will be computed in order to find the target
interval of minimal error penalty values.}
  \item{minutes.limit}{Time limit; the search will stop at a sub-optimal target interval
if this many minutes has elapsed. Useful for testing environments
with build time limits (travis). Default NULL means to use the
value in option PeakSegPipeline.problem.target.minutes (or Inf if
that option is not set).}
  \item{verbose}{
}
}

\value{List of info related to target interval computation: target is the
interval of log(penalty) values that achieve minimum incorrect
labels (numeric vector of length 2), target.iterations is a
data.table with target intervals as a function of iteration,
models is a data.table with one row per model for which the label
error was computed.}

\author{Toby Dylan Hocking}




\examples{

library(PeakSegPipeline)
data(Mono27ac)
## Write the Mono27ac data set to disk.
data.dir <- file.path(
  tempfile(),
  "H3K27ac-H3K4me3_TDHAM_BP",
  "samples",
  "Mono1_H3K27ac",
  "S001YW_NCMLS",
  "problems",
  "chr11:60000-580000")
dir.create(data.dir, recursive=TRUE, showWarnings=FALSE)
write.table(
  Mono27ac$labels, file.path(data.dir, "labels.bed"),
  col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(
  Mono27ac$coverage, file.path(data.dir, "coverage.bedGraph"),
  col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

## Compute target interval. Specifying minutes.limit stops the
## optimization after that number of minutes, resulting in an
## imprecise target interval, but saving time (to avoid NOTE on
## CRAN).
target.list <- problem.target(data.dir, minutes.limit=0.05)

## These are all the models computed in order to find the target
## interval.
print(target.list$models[, list(
  penalty, log.penalty=log(penalty), peaks, total.cost, fn, fp, errors)])

## This is the target interval in log(penalty) values.
print(target.list$target)

}
