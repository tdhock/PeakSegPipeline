library(testthat)
library(PeakSegPipeline)
library(data.table)
context("demo")
chromInfo.dt <- fread(text="chr1 100
chr9 100
chr10 100")
peaks.dt <- fread(text="chr1 1 2 3.4
chr9 10 20 24.81
chr10 30 40 9.20")
chromInfo.txt <- file.path(tempdir(), "chromInfo.txt")
peaks.bedGraph <- file.path(tempdir(), "peaks.bedGraph")
peaks.bigWig <- file.path(tempdir(), "peaks.bigWig")
fwrite(chromInfo.dt, chromInfo.txt, sep="\t", col.names=FALSE)
fwrite(peaks.dt, peaks.bedGraph, sep="\t", col.names=FALSE)

test_that("bedGraphToBigWig sorts alphabetically", {
  bigWig.created <- bedGraphToBigWig(
    peaks.bedGraph, chromInfo.txt, peaks.bigWig)
  expect_identical(bigWig.created, TRUE)
})
