library(testthat)
context("convert labels")
library(PeakSegPipeline)

labels.content <- "
chr2L:72,965-73,182 peakStart
chr2L:73,183-74,241 peaks
chr2L:74,242-74,305 peakEnd

chr2L:11,011,984-11,012,154 peakStart
chr2L:11,012,278-11,012,441 peakEnd
"
proj.dir <- tempfile()
labels.txt <- file.path(proj.dir, "labels", "someLabels.txt")
dir.create(dirname(labels.txt), recursive=TRUE)
writeLines(labels.content, labels.txt)

test_that("informative error for convert_labels on labels file", {
  expect_error({
    convert_labels(labels.txt)
  }, "should be a project directory with proj.dir/labels/*.txt files", fixed=TRUE)
})

test_that("informative error for positive label with no sample groups", {
  expect_error({
    convert_labels(proj.dir)
  }, "need at least one sample group up for each positive label")
})
