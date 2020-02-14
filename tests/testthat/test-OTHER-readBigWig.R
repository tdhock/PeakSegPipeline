library(testthat)
context("readBigWig")

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
coverage.bigWig <- tempfile()
PeakSegPipeline::bedGraphToBigWig(coverage.bedGraph, sizes.tsv, coverage.bigWig)
cov.dt <- sizes.dt[, PeakSegPipeline:::readBigWig(
  coverage.bigWig, chrom, 0, bases)]
test_that("readBigWig returns integer count", {
  expect_is(cov.dt$count, "integer")
})
