library(testthat)
library(PeakSegPipeline)
library(data.table)
context("input")
options(
  mc.cores=parallel::detectCores(),
  PeakSegPipeline.problem.target.minutes=1)

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

## Download bigWig files from github.
bigWig.part.vec <- c(
  "Input/MS010302",
  ## "bcell/MS010302",
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
chr10:35,261,418-35,314,654 peakStart kidney
chr10:35,343,031-35,398,459 peakEnd kidney

chr10:38,041,023-38,102,554 noPeaks
chr10:38,296,008-38,307,179 peakStart kidney
chr10:38,379,045-38,391,967 peakStart kidney
chr10:38,404,899-38,412,089 peakEnd kidney
chr10:38,413,073-38,444,133 noPeaks

chr10:38,585,584-38,643,190 noPeaks
chr10:38,643,191-38,650,766 peakStart kidney
chr10:38,731,066-38,750,574 peakEnd kidney
chr10:38,750,960-38,790,663 noPeaks

chr10:38,807,475-38,815,200 noPeaks
chr10:38,815,201-38,816,355 peakStart kidney Input
chr10:38,818,377-38,819,342 peakEnd kidney Input
"

test.data.dir <- file.path(Sys.getenv("HOME"), "PeakSegPipeline-test")
##test.data.dir <- file.path(tempdir(), "PeakSegPipeline-test")
non.integer.dir <- file.path(test.data.dir, "non-integer")
demo.dir <- file.path(test.data.dir, "input")
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
  demo.bigWig <- sub("non-integer", "input", bigWig.file)
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
  cat("chr10	18024675	38818835", file=problems.bed)
}

## Pipeline should raise error for non-integer data.
system(paste("bigWigToBedGraph", bigWig.file, "/dev/stdout|head"))
test_that("error for non-integer data in bigWigs", {
  expect_error({
    pipeline(non.integer.dir)
  }, "non-integer data in")
})
unlink(non.integer.dir, recursive=TRUE, force=TRUE)

## Pipeline should run to completion for integer count data.
system(paste("bigWigToBedGraph", demo.bigWig, "/dev/stdout|head"))
index.html <- file.path(demo.dir, "index.html")
unlink(index.html)
pipeline(demo.dir)
test_that("index.html is created", {
  expect_true(file.exists(index.html))
})

