library(testthat)
library(data.table)
chrom <- "chr1"
d <- function(chromStart, chromEnd, count){
  data.frame(
    chrom,
    chromStart=as.integer(chromStart),
    chromEnd=as.integer(chromEnd),
    count=as.integer(count))
}
dup.gap <- rbind(
  d(0, 10, 5),
  d(10, 15, 5),
  d(20, 30, 2))
f <- tempfile()
f.bedGraph <- paste0(f, ".bedGraph")
f.bigWig <- paste0(f, ".bigWig")
f.chromInfo <- paste0(f, ".chromInfo")
data.table::fwrite(dup.gap, f.bedGraph, sep="\t", col.names=FALSE)
last <- 50
data.table::fwrite(
  data.frame(chrom, last), f.chromInfo, sep="\t", col.names=FALSE)
PeakSegPipeline::bedGraphToBigWig(f.bedGraph, f.chromInfo, f.bigWig)
f.new <- paste0(f, ".new")
test_that("NoGaps removes gaps and duplicates", {
  PeakSegPipeline::bigWigToBedGraphNoGaps(f.bigWig, chrom, 0, last, f.new)
  cleaned <- data.table::fread(f.new, col.names = names(dup.gap))
  expected <- data.table(rbind(
    d(0, 15, 5),
    d(15, 20, 0),
    d(20, 30, 2),
    d(30, 50, 0)))
  expect_identical(cleaned, expected)
})
