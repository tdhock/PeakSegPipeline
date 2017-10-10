library(testthat)
context("problem.PeakSegFPOP")
library(PeakSegPipeline)

data(chrXcrash)

coverage <- data.frame(chrom="chrX", chrXcrash)
dir.create(prob.dir <- tempfile())
coverage.bedGraph <- file.path(prob.dir, "coverage.bedGraph")
fwrite(coverage, coverage.bedGraph, col.names=FALSE, sep="\t")

fit <- problem.PeakSegFPOP(prob.dir, "0")
