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
test_that("informative error for no sample directories", {
  expect_error({
    convert_labels(proj.dir)
  }, "no sample directories")
})
dir.create(file.path(proj.dir, "samples", "group1", "sample1"), recursive=TRUE)
test_that("informative error for positive label with no sample groups", {
  expect_error({
    convert_labels(proj.dir)
  }, "need at least one sample group up for each positive label")
})

labels.content <- "
chr2L:72,965-73,182 peakStart monocyte
chr2L:73,183-74,241 peaks monocyte
chr2L:74,242-74,305 peakend monocyte

chr2L:11,011,984-11,012,154 peakStart bcell
chr2L:11,012,278-11,012,441 peakEnd bcell
"
proj.dir <- tempfile()
labels.txt <- file.path(proj.dir, "labels", "someLabels.txt")
dir.create(dirname(labels.txt), recursive=TRUE)
dir.create(file.path(proj.dir, "samples", "group1", "sample1"), recursive=TRUE)
writeLines(labels.content, labels.txt)
test_that("informative error for unrecognized annotation", {
  expect_error({
    convert_labels(proj.dir)
  }, "peakEnd")
})

labels.content <- "
chr2L:72,965-73,182 peakStart monocyte
chr2L:73,183-74,241 peaks monocyte
chr2L:74,242-74,305 peakend monocyte

chr3L:9304437-9304498 peakStart sjl_kc
chr3L:9304499-9304558 peaks sjl_kc 
chr3L;9304559-9304571 peakEnd sjl_kc

chr2L:11,011,984-11,012,154 peakStart bcell
chr2L:11,012,278-11,012,441 peakEnd bcell

"
proj.dir <- tempfile()
labels.txt <- file.path(proj.dir, "labels", "someLabels.txt")
dir.create(dirname(labels.txt), recursive=TRUE)
dir.create(file.path(proj.dir, "samples", "group1", "sample1"), recursive=TRUE)
writeLines(labels.content, labels.txt)
test_that("informative error for unrecognized annotation", {
  expect_error({
    convert_labels(proj.dir)
  }, "did not match regex")
})
