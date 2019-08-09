library(testthat)
context("coverage")
library(PeakSegPipeline)
library(data.table)
count.dt <- fread("
chrom	chromStart	chromEnd	coverage
chr1	0	10	1
chr1	10	20	2
chr2	30	40	1
")
chromInfo <- data.table(chrom=c("chr1", "chr2"), chromEnd=50)
work.dir <- tempdir()
chromInfo.txt <- file.path(work.dir, "chromInfo.txt")
fwrite(chromInfo, chromInfo.txt, sep="\t", col.names=FALSE)
input.bedGraph <- file.path(work.dir, "input.bedGraph")
fwrite(count.dt, input.bedGraph, col.names=FALSE, sep="\t")
bedGraph.dt <- bedGraphCoverage(input.bedGraph)

test_that("bedGraphCoverage says total coverage is 40", {
  expect_equal(bedGraph.dt$total.coverage, 40)
})

