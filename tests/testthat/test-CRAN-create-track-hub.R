library(testthat)
library(PeakSegPipeline)
library(data.table)
context("create track hub")

test.data.dir <- file.path(Sys.getenv("HOME"), "PeakSegPipeline-test")
non.integer.dir <- file.path(test.data.dir, "non-integer")
demo.dir <- file.path(test.data.dir, "input")
hub.txt <- file.path(demo.dir, "hub.txt")
trackDb.txt <- file.path(demo.dir, "trackDb.txt")

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
  "kidney/MS002201",
  "bcell/MS010302"
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

# Step 0
convert_labels(demo.dir)
# Step 1
problem.dir.glob <- file.path(demo.dir, "samples", "*", "*", "*", "*")
problem.dir.vec <- Sys.glob(problem.dir.glob)
for(problem.dir in problem.dir.vec){
  problem.target(problem.dir)
}
# Step 2
problem.train(demo.dir)
# Step 3
prob.dir <- file.path(demo.dir, "problems", "chr10:18024675-38818835")
problem.pred.cluster.targets(prob.dir)
# Step 4
problem.joint.train(demo.dir)
# Step 5
job <- file.path(demo.dir, "jobs", "1")
problem.joint.predict.job(job)
# Step 6
plot_all(demo.dir)

url <- paste0("http://CHANGE.THIS/~URL/", basename(demo.dir))
email <- "toby.hocking@r-project.org"
genome <- "hg19"
create_track_hub(demo.dir, url, genome, email)

test_that("hub.txt file is created", {
  expect_true(file.exists(hub.txt))
})

test_that("hub.txt file has correct email", {
  hub.vec <- readLines(hub.txt)
  index <- length(hub.vec)
  hub.txt.email <- substring(hub.vec[index], 7)
  expect_equal(hub.txt.email, email)
})

test_that("trackDb.txt file has the correct link to the bigwig files", {
  trackDb.vec <- readLines(trackDb.txt)
  kidney.coverage.line <- grep("samples/kidney/MS002201/coverage.bigWig", trackDb.vec)
  kidney.jointpeak.line <- grep("samples/kidney/MS002201/joint_peaks.bigWig", trackDb.vec)
  bcell.coverage.line <- grep("samples/bcell/MS010302/coverage.bigWig", trackDb.vec)
  bcell.jointpeak.line <- grep("samples/bcell/MS010302/joint_peaks.bigWig", trackDb.vec)
  split.index <- 14
  trackDb.kidney.coverage.url <- substring(trackDb.vec[kidney.coverage.line], split.index)
  trackDb.kidney.jointpeak.url <- substring(trackDb.vec[kidney.jointpeak.line], split.index)
  trackDb.bcell.coverage.url <- substring(trackDb.vec[bcell.coverage.line], split.index)
  trackDb.bcell.jointpeak.url <- substring(trackDb.vec[bcell.jointpeak.line], split.index)
  sorted.bigWig.part.vec <- sort(bigWig.part.vec)
  bcell.url <- paste0(url, "/samples/", sorted.bigWig.part.vec[1], "/")
  kidney.url <- paste0(url, "/samples/", sorted.bigWig.part.vec[2], "/")
  expect_equal(trackDb.kidney.coverage.url, paste0(kidney.url, "coverage.bigWig"))
  expect_equal(trackDb.kidney.jointpeak.url, paste0(kidney.url, "joint_peaks.bigWig"))
  expect_equal(trackDb.bcell.coverage.url, paste0(bcell.url, "coverage.bigWig"))
  expect_equal(trackDb.bcell.jointpeak.url, paste0(bcell.url, "joint_peaks.bigWig"))
})
