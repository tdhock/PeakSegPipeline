library(testthat)
context("readBigWig")

## This data set has some really large integer values > 1,000,000
## which bigWigToBedGraph prints as 1.0001e6 so fread returns a
## numeric column, which would cause an error if passed to
## PeakSegJoint::ProfileList. (fixed by coercing to integer in
## PeakSegPipeline::problem.joint after fread but before ProfileList)
norm.bedGraph <- system.file(
  "extdata", "norm-big-coverage.bedGraph",
  package="PeakSegPipeline", mustWork=TRUE)
coverage.bedGraph <- tempfile()
PeakSegPipeline:::denormalizeBedGraph(
  norm.bedGraph, coverage.bedGraph)
coverage.dt <- data.table::fread(
  coverage.bedGraph,
  col.names=c("chrom", "chromStart", "chromEnd", "count"))
sizes.dt <- coverage.dt[, .(bases=max(chromEnd)), by=chrom]
sizes.tsv <- tempfile()
data.table::fwrite(sizes.dt, sizes.tsv, sep="\t", col.names=FALSE)

data.dir <- tempfile()
sample.dir <- file.path(data.dir, "samples", "group1", "sample1")
dir.create(sample.dir, showWarnings=FALSE, recursive=TRUE)
coverage.bigWig <- file.path(sample.dir, "coverage.bigWig")
PeakSegPipeline::bedGraphToBigWig(coverage.bedGraph, sizes.tsv, coverage.bigWig)
prob.dt <- coverage.dt[, .(
  chrom=chrom[1],
  start=chromStart[1],
  end=chromEnd[.N])]
prob.name <- prob.dt[, paste0(chrom, ":", start, "-", end)]
jprob.dir <- file.path(
  data.dir, "problems", prob.name, "jointProblems", prob.name)
dir.create(jprob.dir, showWarnings=FALSE, recursive=TRUE)
test_that("problem.joint returns list", {
  fit <- PeakSegPipeline:::problem.joint(jprob.dir)
  expect_is(fit, "list")
})
