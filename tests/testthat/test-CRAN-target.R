library(testthat)
context("problem.target")
library(PeakSegPipeline)

data(Mono27ac, envir=environment())
## Write the Mono27ac data set to disk.
problem.dir <- file.path(
  tempfile(),
  "H3K27ac-H3K4me3_TDHAM_BP",
  "samples",
  "Mono1_H3K27ac",
  "S001YW_NCMLS",
  "problems",
  "chr11:60000-580000")
dir.create(problem.dir, recursive=TRUE, showWarnings=FALSE)
problem.labels.bed <- file.path(problem.dir, "labels.bed")
write.table(
  Mono27ac$labels, problem.labels.bed,
  col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(
  Mono27ac$coverage, file.path(problem.dir, "coverage.bedGraph"),
  col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(
  data.frame(minutes=0.05), file.path(problem.dir, "target.minutes"),
  col.names=FALSE, row.names=FALSE, quote=FALSE)
test_that("problem.target uses problem/labels.bed if present", {
  target.list <- problem.target(problem.dir)
  expect_is(target.list, "list")
  expect_is(target.list$target, "numeric")
  expect_identical(length(target.list$target), 2L)
  expect_is(target.list$target.iterations, "data.table")
  expect_is(target.list$models, "data.table")
})

unlink(problem.labels.bed)
test_that("problem.target errors if no labels", {
  expect_error({
    problem.target(problem.dir)
  }, "need labels to compute target interval")
})

problems.dir <- dirname(problem.dir)
sample.dir <- dirname(problems.dir)
sample.labels.bed <- file.path(problem.dir, "labels.bed")
write.table(
  Mono27ac$labels, sample.labels.bed,
  col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
test_that("problem.target uses sampleID/labels.bed if present", {
  target.list <- problem.target(problem.dir)
  expect_is(target.list, "list")
  expect_is(target.list$target, "numeric")
  expect_identical(length(target.list$target), 2L)
  expect_is(target.list$target.iterations, "data.table")
  expect_is(target.list$models, "data.table")
})
