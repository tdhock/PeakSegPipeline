library(testthat)
library(PeakSegPipeline)
context("gap2problems")
library(data.table)

gaps <- data.table(
  chrom="chr1",
  gapStart=c(10, 30, 40),
  gapEnd=c(20, 40, 50))
fwrite(gaps, gaps.txt <- file.path(tempdir(), "gaps.txt"))
sizes <- data.table(
  chrom=c("chr1", "chrM"),
  bases=c(200, 100))
fwrite(sizes, sizes.txt <- file.path(tempdir(), "sizes.txt"))
gap2problems(gaps.txt, sizes.txt, contigs.txt <- file.path(tempdir(), "contigs.txt"))
contigs <- fread(contigs.txt, col.names=c("chrom", "problemStart", "problemEnd"))[order(chrom)]

test_that("gap2problems creates contig for chrM even though there are no gaps", {
  expect_identical(contigs$chrom, c("chr1", "chr1", "chr1", "chrM"))
  expect_identical(contigs$problemStart, as.integer(c(0, 20, 50, 0)))
  expect_identical(contigs$problemEnd, as.integer(c(10, 30, 200, 100)))
})
