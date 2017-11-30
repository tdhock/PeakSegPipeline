library(testthat)
context("denormalize")
library(PeakSegPipeline)
library(data.table)
norm.dt <- fread("
chrom	chromStart	chromEnd	coverage
chr1	9979	9993	0.07
chr1	9993	9998	0.14
chr1	9998	10024	0.21
chr1	10024	10043	0.28
chr1	10043	10049	0.35
chr1	10049	10054	0.42
chr1	10054	10091	0.48
chr1	10091	10135	0.55
chr1	10135	10149	0.48
chr1	10149	10154	0.42
")
chromInfo <- fread("
chr1 100000
")
work.dir <- tempdir()
chromInfo.txt <- file.path(work.dir, "chromInfo.txt")
fwrite(chromInfo, chromInfo.txt, sep="\t", col.names=FALSE)
input.bedGraph <- file.path(work.dir, "input.bedGraph")
fwrite(norm.dt, input.bedGraph, col.names=FALSE, sep="\t")
norm.bigWig <- file.path(work.dir, "norm.bigWig")
system.or.stop(paste(
  "bedGraphToBigWig",
  input.bedGraph,
  chromInfo.txt,
  norm.bigWig))

denorm.bigWig <- file.path(work.dir, "denorm.bigWig")
denormalizeBigWig(norm.bigWig, denorm.bigWig)

test_that("intermediate bedGraph files are deleted", {
  denorm.bedGraph <- file.path(work.dir, "denorm.bedGraph")
  norm.bedGraph <- file.path(work.dir, "norm.bedGraph")
  expect_true(!file.exists(denorm.bedGraph))
  expect_true(!file.exists(norm.bedGraph))
})

denorm.dt <- readBigWig(denorm.bigWig, "chr1", 0, 100000)
test_that("coverage converted to integers", {
  expect_true(is.integer(denorm.dt$count))
})

## second test, 0.1 0.2 etc.
norm.dt <- fread("
chrom	chromStart	chromEnd	coverage
chr1	9993	9998	0.1
chr1	9998	10024	0.2
chr1	10024	10043	0.3
chr1	10043	10049	0.4
chr1	10049	10054	0.5
chr1	10054	10091	0.6
chr1	10091	10135	0.7
chr1	10135	10149	0.8
chr1	10149	10154	0.9
")
chromInfo <- fread("
chr1 100000
")
work.dir <- tempdir()
chromInfo.txt <- file.path(work.dir, "chromInfo.txt")
fwrite(chromInfo, chromInfo.txt, sep="\t", col.names=FALSE)
input.bedGraph <- file.path(work.dir, "input.bedGraph")
fwrite(norm.dt, input.bedGraph, col.names=FALSE, sep="\t")
norm.bigWig <- file.path(work.dir, "norm.bigWig")
system.or.stop(paste(
  "bedGraphToBigWig",
  input.bedGraph,
  chromInfo.txt,
  norm.bigWig))
denorm.bigWig <- file.path(work.dir, "denorm.bigWig")
denormalizeBigWig(norm.bigWig, denorm.bigWig)
denorm.dt <- readBigWig(denorm.bigWig, "chr1", 0, 100000)
test_that("coverage counts 0:9", {
  expect_identical(denorm.dt$count, 1:9)
})
