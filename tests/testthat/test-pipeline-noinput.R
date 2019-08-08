library(testthat)
library(PeakSegPipeline)
library(data.table)
context("noinput")
options(
  mc.cores=parallel::detectCores())
test.data.dir <- file.path(Sys.getenv("HOME"), "PeakSegPipeline-test")
non.integer.dir <- file.path(test.data.dir, "non-integer")
demo.dir <- file.path(test.data.dir, "noinput")
index.html <- file.path(demo.dir, "index.html")
download.to <- function
(u, f, writeFun=if(grepl("bigWig", f))writeBin else writeLines){
  if(!file.exists(f)){
    require(httr)
    f.dir <- dirname(f)
    dir.create(f.dir, showWarnings=FALSE, recursive=TRUE)
    request <- GET(u)
    stop_for_status(request)
    writeFun(content(request), f)
  }
}
res.list <- list(
  walltime = 3600, #in minutes
  ncpus=1,
  ntasks=1,
  chunks.as.arrayjobs=TRUE)
## Download bigWig files from github.
bigWig.part.vec <- c(
  "Input/MS010302",
  "bcell/MS010302",
  ## "Input/MS002202",
  ## "kidney/MS002202",
  ## "Input/MS026601",
  ## "bcell/MS026601",
  ## "Input/MS002201",
  "kidney/MS002201"
)

label.txt <- "
chr10:33,061,897-33,162,814 noPeaks
chr10:33,456,000-33,484,755 peakStart kidney
chr10:33,597,317-33,635,209 peakEnd kidney
chr10:33,662,034-33,974,942 noPeaks

chr10:35,182,820-35,261,001 noPeaks
chr10:35,261,418-35,314,654 peakStart bcell kidney
chr10:35,343,031-35,398,459 peakEnd bcell kidney

chr10:38,041,023-38,102,554 noPeaks
chr10:38,296,008-38,307,179 peakStart bcell kidney
chr10:38,379,045-38,391,967 peakStart bcell kidney
chr10:38,404,899-38,412,089 peakEnd bcell kidney
chr10:38,413,073-38,444,133 noPeaks

chr10:38,585,584-38,643,190 noPeaks
chr10:38,643,191-38,650,766 peakStart bcell kidney
chr10:38,731,066-38,750,574 peakEnd bcell kidney
chr10:38,750,960-38,790,663 noPeaks
"

chrom.sizes.file <- tempfile()
chrom.sizes <- data.table(chrom="chr10", bases=128616069)
fwrite(chrom.sizes, chrom.sizes.file, sep="\t", col.names=FALSE)
repos.url <- "https://raw.githubusercontent.com/tdhock/input-test-data/master/"
for(bigWig.part in bigWig.part.vec){
  bigWig.file <- file.path(
    non.integer.dir, "samples",
    bigWig.part, "coverage.bigWig")
  bigWig.url <- paste0(repos.url, bigWig.part, ".bigwig")
  download.to(bigWig.url, bigWig.file)
  demo.bigWig <- sub("non-integer", "noinput", bigWig.file)
  if(!file.exists(demo.bigWig)){
    dir.create(dirname(demo.bigWig), showWarnings=FALSE, recursive=TRUE)
    bw.dt <- readBigWig(bigWig.file, "chr10", 0, 128616069)
    out.dt <- data.table(chrom="chr10", bw.dt)
    demo.bedGraph <- sub("bigWig", "bedGraph", demo.bigWig)
    fwrite(out.dt, demo.bedGraph, sep="\t", col.names=FALSE)
    system.or.stop(
      paste("bedGraphToBigWig", demo.bedGraph, chrom.sizes.file, demo.bigWig))
    unlink(demo.bedGraph)
  }
}

for(set.dir in c(non.integer.dir, demo.dir)){
  labels.file <- file.path(set.dir, "labels", "some_labels.txt")
  dir.create(dirname(labels.file), showWarnings=FALSE, recursive=TRUE)
  writeLines(label.txt, labels.file)
  problems.bed <- file.path(set.dir, "problems.bed")
  unlink(problems.bed)
  cat("chr10	60000	17974675
chr10	18024675	38818835
chr10	38868835	39154935
chr10	42746000	46426964
chr10	47529169	47792476
chr10	47892476	48055707
chr10	48105707	49095536
chr10	49195536	51137410
chr10	51187410	51398845
chr10	51448845	125869472
chr10	125919472	128616069
", file=problems.bed)
}

## Pipeline should raise error for non-integer data.
test_that("error for non-integer data in bigWigs", {
  expect_error({
    pipeline(non.integer.dir)
  }, "non-integer data in")
})
unlink(non.integer.dir, recursive=TRUE, force=TRUE)

## Set time limit.
(sample.dir.vec <- Sys.glob(file.path(
  demo.dir, "samples", "*", "*")))
prob.dir.vec <- file.path(
  sample.dir.vec, "problems", "chr10:18024675-38818835")
limit.dt <- data.table(minutes=5)
for(prob.dir in prob.dir.vec){
  dir.create(prob.dir, showWarnings=FALSE, recursive=TRUE)
  limit.file <- file.path(prob.dir, "target.minutes")
  fwrite(limit.dt, limit.file, col.names=FALSE)
}

## Remove one sampleID/problems dir to simulate what happens when
## running jobs_create (which does not create problems dirs) then
## jobs_submit.
one.problems.dir <- dirname(prob.dir)
unlink(one.problems.dir, recursive=TRUE)
test_that("problem.coverage makes a directory", {
  prob <- problem.coverage(prob.dir)
  expect_true(file.exists(file.path(prob.dir, "coverage.bedGraph")))
})
limit.file <- file.path(prob.dir, "target.minutes")
fwrite(limit.dt, limit.file, col.names=FALSE)

jobs <- jobs_create(demo.dir, verbose=1)
test_that("jobs_create returns dt", {
  expect_identical(names(jobs), c("step", "fun", "arg"))
  expect_is(jobs, "data.table")
})

## Pipeline should run to completion using SLURM.
unlink(index.html)
test_that("index.html is created via batchtools", {
  reg.list <- jobs_submit_batchtools(jobs, res.list)
  reg <- reg.list[[length(reg.list)]]
  result <- batchtools::waitForJobs(reg=reg, sleep=function(i){
    system("squeue")
    10
  })
  expect_true(file.exists(index.html))
  log.glob <- file.path(demo.dir, "registry", "*", "logs", "*")
  system(paste("tail -n 10000", log.glob))
})

test_that("entries of peaks matrix are 0/1", {
  mat.tsv.gz <- file.path(demo.dir, "peaks_matrix_sample.tsv.gz")
  peak.dt <- fread(cmd=paste("zcat", mat.tsv.gz))
  class.vec <- as.character(sapply(peak.dt, class))
  expected.class.vec <- c("character", rep("integer", length(bigWig.part.vec)))
  expect_identical(class.vec, expected.class.vec)
  binary.mat <- as.matrix(peak.dt[, -1, with=FALSE])
  expect_true(all(binary.mat %in% c(0,1)))
})

some.jobs <- jobs[step >= 5]
unlink(index.html)
test_that("run only steps 5-6 creates index.html", {
  reg.list <- jobs_submit_batchtools(some.jobs, res.list)
  reg <- reg.list[[length(reg.list)]]
  result <- batchtools::waitForJobs(reg=reg, sleep=function(i){
    system("squeue")
    10
  })
  expect_true(file.exists(index.html))
  for(step.i in names(reg.list)){
    log.glob <- file.path(demo.dir, "registry", step.i, "logs", "*")
    system(paste("tail -n 10000", log.glob))
  }
})
