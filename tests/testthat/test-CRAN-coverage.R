library(testthat)
context("coverage")
library(PeakSegPipeline)
library(data.table)
count.dt <- fread("
chrom	chromStart	chromEnd	coverage
chr1	0	10	1
chr1	10	20	2
chr1	30	40	1
")
chromInfo <- fread("
chr1 100
")
work.dir <- tempdir()
chromInfo.txt <- file.path(work.dir, "chromInfo.txt")
fwrite(chromInfo, chromInfo.txt, sep="\t", col.names=FALSE)
input.bedGraph <- file.path(work.dir, "input.bedGraph")
fwrite(count.dt, input.bedGraph, col.names=FALSE, sep="\t")
bedGraph.dt <- bedGraphCoverage(input.bedGraph)

test_that("bedGraphCoverage says total coverage is 40", {
  expect_equal(bedGraph.dt$total.coverage, 40)
})

count.bigWig <- file.path(work.dir, "count.bigWig")
system.or.stop(paste(
  "bedGraphToBigWig",
  input.bedGraph,
  chromInfo.txt,
  count.bigWig))

bigWig.dt <- bigWigCoverage(count.bigWig)

test_that("bigWigCoverage stats are OK", {
  expect_equal(bigWig.dt$total.coverage, 40)
  expect_equal(bigWig.dt$mean.coverage, 0.4)
  expect_equal(bigWig.dt$total.bases, 100)
})

test_that("intermediate file is deleted", {
  count.bedGraph <- file.path(work.dir, "count.bedGraph")
  expect_true(!file.exists(count.bedGraph))
})
