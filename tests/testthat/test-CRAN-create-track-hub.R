library(testthat)
library(PeakSegPipeline)
library(data.table)
context("create track hub")

test.data.dir <- file.path(Sys.getenv("HOME"), "PeakSegPipeline-test")
non.integer.dir <- file.path(test.data.dir, "non-integer")
demo.dir <- file.path(test.data.dir, "input")
hub.txt <- file.path(demo.dir,"hub.txt")

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

### Download bigWig files from github.
bigWig.part.vec <- c(
  "Input/MS010302",
  "bcell/MS010302",
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

chr10:38,807,475-38,815,200 noPeaks
chr10:38,815,201-38,816,355 peakStart bcell kidney Input
chr10:38,818,377-38,819,342 peakEnd bcell kidney Input

chr10:39,098,319-39,111,384 noPeaks
chr10:39,125,134-39,125,550 peakStart bcell kidney Input
chr10:39,125,594-39,126,266 peakEnd bcell kidney Input
chr10:39,126,866-39,140,858 noPeaks
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


sample.dir <- dirname(demo.bigWig)
problem.dir <- file.path(sample.dir, "problems", "chr10:18024675-38818835")
coverage.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
unlink(coverage.bedGraph)

for(set.dir in c(non.integer.dir, demo.dir)){
  labels.file <- file.path(set.dir, "labels", "some_labels.txt")
  dir.create(dirname(labels.file), showWarnings=FALSE, recursive=TRUE)
  writeLines(label.txt, labels.file)
  problems.bed <- file.path(set.dir, "problems.bed")
  unlink(problems.bed)
  cat("chr10	18024675	38818835", file=problems.bed)
}

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

system(paste("bigWigToBedGraph", demo.bigWig, "/dev/stdout|head"))

jobs_create_run(demo.dir)

test_that("hub.txt file is created", {
  expect_true(file.exists(hub.txt))
})
