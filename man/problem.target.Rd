\name{problem.target}
\alias{problem.target}
\title{problem target}
\description{Compute target interval for a segmentation problem. This function
repeated calls problem.PeakSegFPOP with different penalty values,
until it finds an interval of penalty values with minimal label
error. The calls to PeakSegFPOP are parallelized using mclapply if
you set options(mc.cores).
A time limit in minutes may be specified in a file
problem.dir/target.minutes;
the search will stop at a sub-optimal target interval
if this many minutes has elapsed. Useful for testing environments
with build time limits (travis). }
\usage{problem.target(problem.dir, verbose = 0)}
\arguments{
  \item{problem.dir}{problemID directory in which coverage.bedGraph has already been
computed. If there is a labels.bed file then the number of
incorrect labels will be computed in order to find the target
interval of minimal error penalty values.}
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
data(Mono27ac, envir=environment())
## Write the Mono27ac data set to disk.
problem.dir <- file.path(
  tempfile(),
  "H3K27ac-H3K4me3_TDHAM_BP",
  "samples",
  "Mono1_H3K27ac",
  "S001YW_NCMLS",
  "problems",
  "chr11:60000-580000")
dir.create(problem.dir, recursive=TRUE, showWarnings=FALSE)
write.table(
  Mono27ac$labels, file.path(problem.dir, "labels.bed"),
  col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(
  Mono27ac$coverage, file.path(problem.dir, "coverage.bedGraph"),
  col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
## Creating a target.minutes file stops the optimization after that
## number of minutes, resulting in an imprecise target interval, but
## saving time (to avoid NOTE on CRAN).
write.table(
  data.frame(minutes=0.05), file.path(problem.dir, "target.minutes"),
  col.names=FALSE, row.names=FALSE, quote=FALSE)

## Compute target interval.
target.list <- problem.target(problem.dir)

## These are all the models computed in order to find the target
## interval.
print(target.list$models[, list(
  penalty, log.penalty=log(penalty), peaks, total.loss, fn, fp, errors)])

## This is the target interval in log(penalty) values.
print(target.list$target)

}
