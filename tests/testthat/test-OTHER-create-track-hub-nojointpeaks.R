library(testthat)
library(PeakSegPipeline)
library(data.table)
context("create track hub with no joint peaks")
test.data.dir <- file.path(Sys.getenv("HOME"), "PeakSegPipeline-test")
non.integer.dir <- file.path(test.data.dir, "non-integer")
demo.dir <- file.path(test.data.dir, "input")
hub.txt <- file.path(demo.dir, "hub.txt")
trackDb.txt <- file.path(demo.dir, "trackDb.txt")
### Download bigWig files from github.
bigWig.part.vec <- c(
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
    cmd <- paste("bedGraphToBigWig", demo.bedGraph, chrom.sizes.file, demo.bigWig)
    system.or.stop(cmd)
    unlink(demo.bedGraph)
  }
}
unlink(non.integer.dir, recursive=TRUE, force=TRUE)
url <- paste0("http://CHANGE.THIS/~URL/", basename(demo.dir))
email <- "toby.hocking@r-project.org"
genome <- "hg19"
create_track_hub(demo.dir, url, genome, email)

test_that("hub.txt file is created", {
  expect_true(file.exists(hub.txt))
})

test_that("hub.txt file has correct email", {
  hub.vec <- readLines(hub.txt)
  email.line <- grep(email, hub.vec)
  expected.email <- sub(".* ", "", hub.vec[email.line])
  expect_equal(email, expected.email)
})

test_that("trackDb.txt file has the correct link to the bigwig files", {
  trackDb.vec <- readLines(trackDb.txt)
  coverage.line <- grep("samples/kidney/MS002201/coverage.bigWig", trackDb.vec)
  coverage.url <- sub(".* ", "", trackDb.vec[coverage.line])
  expected.coverage.url <- paste0(url, "/samples/", bigWig.part.vec[1], "/", "coverage.bigWig")
  expect_equal(coverage.url, expected.coverage.url)
})
