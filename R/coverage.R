bigWigCoverage <- function
### Compute total and average coverage from a bigWig file.
(input.bigWig
### Input bigWig file.
){
  stopifnot(is.character(input.bigWig))
  input.path <- normalizePath(input.bedGraph, mustWork=TRUE)
  input.bedGraph <- sub("bigWig$", "bedGraph", input.bigWig)
  if(input.bedGraph==input.bigWig){
    stop("input.bigWig must be a filename that ends with bigWig")
  }
  system.or.stop(paste(
    "bigWigToBedGraph",
    input.bigWig,
    input.bedGraph))
  cov.dt <- bedGraphCoverage(input.bedGraph)
  info.dt <- bigWigInfo(input.bigWig)
  unlink(input.bedGraph)
  cov.dt[, total.bases := sum(info.dt$chromEnd)]
  cov.dt[, mean.coverage := total.coverage/total.bases]
  data.table(input.bigWig, cov.dt)
### one row data.table with columns total.bases, mean.coverage,
### total.coverage, and input.bedGraph.
}

bedGraphCoverage <- function
### Compute total coverage from a bedGraph file.
(input.bedGraph
### Input bedGraph file.
){
  stopifnot(is.character(input.bedGraph))
  input.path <- normalizePath(input.bedGraph, mustWork=TRUE)
  cov.vec <- -1.0
  result <- .C(
    "coverage_interface",
    bedGraph=input.path,
    coverage=cov.vec,
    PACKAGE="PeakSegPipeline")
  data.table(
    input.bedGraph,
    total.coverage=result$coverage)
### data.table of file name and total coverage.
}

