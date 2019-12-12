library(testthat)
library(PeakSegPipeline)
library(data.table)
context("input")
test.data.dir <- file.path(Sys.getenv("HOME"), "PeakSegPipeline-test")
non.integer.dir <- file.path(test.data.dir, "non-integer")
demo.dir <- file.path(test.data.dir, "input")
index.html <- file.path(demo.dir, "index.html")
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
    denormalizeBigWig(bigWig.file, demo.bigWig)
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
test_that("error for non-integer data in bigWigs", {
  expect_error({
    jobs_create_run(non.integer.dir)
  }, "non-integer data in")
})
unlink(non.integer.dir, recursive=TRUE, force=TRUE)

## test for informative error early if ucsc not available.
path.vec <- stop.without.ucsc()
prog <- path.vec[["bigWigInfo"]]
old.mode <- file.info(prog)$mode
Sys.chmod(prog, "444") #read, not write, not exec.
test_that("pipeline fails early if UCSC not available", {
  expect_error({
    jobs_create_run(demo.dir)
  }, "bigWigInfo")
})
Sys.chmod(prog, old.mode)

## Pipeline should run to completion for integer count data.
unlink(index.html)
test_that("index.html is created", {
  jobs_create_run(demo.dir, target.minutes=5)
  expect_true(file.exists(index.html))
})

test_that("relatives links for images", {
  index.vec <- readLines(index.html)
  index.txt <- paste(index.vec, collapse="\n")
  f <- function(x)nc::field(x, '="', '[^"]+')
  match.dt <- nc::capture_all_str(
    index.txt,
    '<a ',
    f("href"),
    "[^<]+",
    '<img ',
    f("src"))
  chunk.limits.csv <- file.path(demo.dir, "chunk.limits.csv")
  chunk.dt <- fread(chunk.limits.csv)
  prefix.vec <- chunk.dt[, paste0(
    "problems/chr10:18024675-38818835/chunks/",
    chrom, ":", chromStart, "-", chromEnd,
    "/")]
  src.vec <- paste0(prefix.vec, "figure-predictions-thumb.png")
  expect_identical(match.dt$src, src.vec)
  href.vec <- paste0(prefix.vec, "figure-predictions.png")
  expect_identical(match.dt$href, href.vec)
})

test_that("joint_peaks.bigWig files have the right number of peaks", {
  jobPeaks.RData.vec <- Sys.glob(file.path(
    demo.dir, "jobs", "*", "jobPeaks.RData"))
  peak.mat.list <- list()
  for(jobPeaks.RData in jobPeaks.RData.vec){
    load(jobPeaks.RData)
    peak.mat.list[[jobPeaks.RData]] <- do.call(cbind, jobPeaks$sample.peaks.vec)
  }
  peak.mat <- do.call(cbind, peak.mat.list)
  library(Matrix)#for rowSums.Matrix
  expected.peaks <- rowSums(peak.mat)
  observed.peaks <- expected.peaks
  for(sample.path in names(expected.peaks)){
    peaks.bigWig <- file.path(
      demo.dir, "samples", sample.path, "joint_peaks.bigWig")
    peaks.dt <- readBigWig(peaks.bigWig, "chr10", 0, 135534747)
    observed.peaks[[sample.path]] <- nrow(peaks.dt)
  }
  expect_equal(observed.peaks, expected.peaks)
})

test_that("no duplicate models in problem cache", {
  models.rds <- Sys.glob(file.path(
    demo.dir, "samples", "*", "*", "problems", "*", "models.rds"))[1]
  prob.dir <- dirname(models.rds)
  models.dt <- PeakSegPipeline:::problem.models(prob.dir)
  count.tab <- table(table(models.dt$penalty.str))
  expect_identical(names(count.tab), "1")
  inf.fit <- PeakSegDisk::PeakSegFPOP_dir(prob.dir, "Inf")
  models.new <- PeakSegPipeline:::problem.models(prob.dir)
  count.new <- table(table(models.new$penalty.str))
  expect_identical(names(count.new), "1")
})

test_that("no duplicate observations in train data cache", {
  train_data.csv <- file.path(demo.dir, "train_data.csv")
  train.orig <- fread(file=train_data.csv)[order(problem.dir)]
  count.orig <- table(table(train.orig$problem.dir))
  expect_identical(names(count.orig), "1")
  tlist <- problem.target(train.orig$problem.dir[1])
  problem.train(demo.dir)
  train.new <- fread(file=train_data.csv)[order(problem.dir)]
  expect_identical(train.new$problem.dir, train.orig$problem.dir)
  count.new <- table(table(train.new$problem.dir))
  expect_identical(names(count.new), "1")
})

test_that("problem.target does not waste time on very similar penalties", {
  problem.dir <- normalizePath(file.path(
    demo.dir, "samples/kidney/MS002201/problems/chr10:18024675-38818835"))
  tlist.1 <- problem.target(problem.dir)
  max.ok <- 2
  ## First remove any models that are duplicated.
  counts.1 <- tlist.1$models[, .(penalties=.N), by=peaks]
  keep.peaks <- counts.1[penalties<=max.ok, peaks]
  keep.models <- tlist.1$models[peaks %in% keep.peaks]
  saveRDS(keep.models, file.path(problem.dir, "models.rds"))
  ## Now run the algo again.
  tlist.2 <- problem.target(problem.dir)
  tlist.2$models[, comp.before := penalty %in% tlist.1$models$penalty]
  print(tlist.2$models[order(penalty), .(penalty, peaks, comp.before)])
  counts.2 <- tlist.2$models[, .(penalties=.N), by=peaks]
  too.many.2 <- counts.2[max.ok<penalties]
  expect_equal(nrow(too.many.2), 0)
  ## Now run again... buggy version got more models here.
  tlist.3 <- problem.target(problem.dir)
  tlist.3$models[, comp.before := penalty %in% tlist.2$models$penalty]
  print(tlist.3$models[order(penalty), .(penalty, peaks, comp.before)])
  counts.3 <- tlist.3$models[, .(penalties=.N), by=peaks]
  too.many.3 <- counts.3[max.ok<penalties]
  expect_equal(nrow(too.many.3), 0)
})
