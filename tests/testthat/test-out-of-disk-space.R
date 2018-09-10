library(testthat)
context("out of disk space")
library(PeakSegPipeline)
library(data.table)
data(Mono27ac)

## sudo mount -t tmpfs -o size=256K tmpfs /tmp/tmp256K

## This test relies on the existence of this very small 256K
## filesystem mounted on /tmp/tmp256K
unlink("/tmp/tmp256K/*", recursive=TRUE)
data.dir <- file.path(
  "/tmp/tmp256K",
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

test_that("problem.PeakSegFPOP errors without crashing R", {
  expect_error({
    problem.PeakSegFPOP(data.dir, "0")
  })
})

