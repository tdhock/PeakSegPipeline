library(testthat)
context("problem.PeakSegFPOP")
library(PeakSegPipeline)

library(data.table)
cov.dt <- fread("chrom   chromStart chromEnd count
chr6_dbb_hap3   3491790 3491834 2
chr6_dbb_hap3   3491834 3491836 1
chr6_dbb_hap3   3491836 3697362 0
chr6_dbb_hap3   3697362 3697408 1
chr6_dbb_hap3   3697408 3701587 0
chr6_dbb_hap3   3701587 3701633 1
chr6_dbb_hap3   3701633 3736386 0
")
prob.dir <- file.path(
  tempfile(),
  "samples",
  "sample1",
  "problems",
  "chr6_dbb_hap3:3491790-3736386")
dir.create(prob.dir, showWarnings=FALSE, recursive=TRUE)
fwrite(
  cov.dt,
  file.path(prob.dir, "coverage.bedGraph"),
  sep="\t", row.names=FALSE, col.names=FALSE)

test_that("no need for problem.bed when running problem.coverage", {
  problem.coverage(prob.dir)
})

test_that("large penalty should not crash solver", {
  fit <- problem.PeakSegFPOP(prob.dir, "866939314852865280")
  expect_identical(fit$loss$peaks, 0L)
})
