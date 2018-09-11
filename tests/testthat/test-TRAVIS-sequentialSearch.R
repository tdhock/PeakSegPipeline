library(testthat)
context("problem.sequentialSearch")
library(PeakSegPipeline)
library(data.table)
data(Mono27ac)

data.dir <- file.path(
  tempfile(),
  "H3K27ac-H3K4me3_TDHAM_BP",
  "samples",
  "Mono1_H3K27ac",
  "S001YW_NCMLS",
  "problems",
  "chr11:60000-580000")
dir.create(data.dir, recursive=TRUE, showWarnings=FALSE)
labels.bed <- file.path(data.dir, "labels.bed")
coverage.bedGraph <- file.path(data.dir, "coverage.bedGraph")
fwrite(
  Mono27ac$labels, labels.bed,
  col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
fwrite(
  Mono27ac$coverage, coverage.bedGraph,
  col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
fwrite(
  data.table(chrom="chr11", chromStart=60000, chromEnd=580000),
  file.path(data.dir, "problem.bed"),
  col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

fit <- problem.sequentialSearch(data.dir, 19L, verbose=1)
test_that("sequential search finds 19 peaks", {
  expect_is(fit$others, "data.table")
  expect_identical(fit$loss$peaks, 19L)
})

most.peaks <- fit$others[penalty==0, peaks]
most.fit <- problem.sequentialSearch(data.dir, most.peaks-1L, verbose=1)
test_that("sequential search returns something for max_peaks-1", {
  expect_is(most.fit$others, "data.table")
})
