library(testthat)
context("out of disk space")
library(PeakSegPipeline)
library(data.table)
data(Mono27ac)

## sudo mount -t tmpfs -o size=256K tmpfs /tmp/tmp256K

## This test relies on the existence of this very small 256K
## filesystem mounted on /tmp/tmp256K
tmp.dir <- "/tmp/tmp256K"
unlink(file.path(tmp.dir, "*"), recursive=TRUE)
data.dir <- file.path(
  tmp.dir,
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
du.dt <- fread(paste("du -bs", file.path(data.dir, "*")))
sum(du.dt$V1)
system(paste("du -bs", tmp.dir))

test_that("problem.PeakSegFPOP error writing cost function database", {
  expect_error({
    problem.PeakSegFPOP(data.dir, "0")
  }, "unable to write to cost function database file")
})

## sudo mount -t tmpfs -o size=4300K tmpfs /tmp/tmp4300K
size <- "160000"
tmp.dir <- paste0("/tmp/tmp", size)
cmd <- paste0(
  "mkdir ",
  tmp.dir,
  "&& sudo mount -t tmpfs -o size=",
  size,
  " tmpfs ",
  tmp.dir)
## system(cmd)
## This test should gracefully fail writing the segments file during
## decoding.
unlink(file.path(tmp.dir, "*"), recursive=TRUE)
data.dir <- file.path(
  tmp.dir,
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

test_that("problem.PeakSegFPOP error writing loss output", {
  expect_error({
    L <- PeakSegFPOP_disk(file.path(data.dir, "coverage.bedGraph"), "Inf")
  }, "unable to write to loss output file")
})

du.dt <- fread(paste("du -bs", file.path(data.dir, "*")))
setnames(du.dt, c("bytes", "file"))
du.dt[, list(bytes, file=sub(".*/", "", file))]
