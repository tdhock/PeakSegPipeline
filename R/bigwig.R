### Run fread but do not stop for an error on an empty file.
fread.or.null <- function(...){
  tryCatch({
    fread(...)
  }, error=function(e){
    NULL
  })
}

readBigWig <- function
### Read part of a bigWig file into R as a data.table (assumes
### bigWigToBedGraph is present on your PATH).
(bigwig.file,
### path or URL of bigwig file.
 chrom,
### chromosome to read.
 start,
### position before reading.
 end
### plain text file where coverage is saved before reading into R.
 ){
  count <- . <- chromStart <- chromEnd <- NULL
  stopifnot(length(bigwig.file) == 1)
  stopifnot(length(chrom) == 1)
  stopifnot(length(start) == 1)
  stopifnot(length(end) == 1)
  stopifnot(is.character(bigwig.file))
  stopifnot(is.character(chrom))
  start <- as.integer(start)
  end <- as.integer(end)
  stopifnot(0 <= start)
  stopifnot(start < end)
  stopifnot(end < Inf)
  cmd <-
    sprintf("bigWigToBedGraph -chrom=%s -start=%d -end=%d %s /dev/stdout",
            chrom, start, end,
            bigwig.file)
  bg <- fread.or.null(cmd, drop=1)
  if(is.null(bg)){
    data.table(chromStart=integer(),
               chromEnd=integer(),
               count=integer())
  }else{
    setnames(bg, c("chromStart", "chromEnd", "norm"))
    stopifnot(0 <= bg$norm)
    nonzero <- bg[0 < norm, ]
    min.nonzero.norm <- min(nonzero[, norm])
    nonzero[, count := as.integer(norm/min.nonzero.norm) ]
    nonzero[, .(
      chromStart,
      chromEnd,
      count
      )]
  }
### data.table with columns chromStart chromEnd count.
}

bigWigInfo <- function
### Run bigWigInfo to find chrom sizes.
(bigwig.file
### path or URL of bigwig file.
 ){
  stopifnot(is.character(bigwig.file))
  stopifnot(length(bigwig.file) == 1)
  cmd <- paste("bigWigInfo", bigwig.file, "-chroms | grep '^\\s'")
  chroms <- fread(cmd, header=FALSE, sep=" ")
  setnames(chroms, c("chrom", "chromStart", "chromEnd"))
  chroms$chrom <- sub("\\s*", "", chroms$chrom)
  chroms
}

### read a bunch of bigwig files into R as a list of data.frames that
### can be passed to PeakSegJointSeveral.
readBigWigSamples <- function(problem, bigwig.file.vec){
  counts.by.sample <- list()
  for(sample.id in names(bigwig.file.vec)){
    bigwig.file <- bigwig.file.vec[[sample.id]]
    sample.counts <- with(problem, {
      readBigWig(bigwig.file, chrom,
                 problemStart, problemEnd)
    })
    ## Make a data.frame and not a data.table, since we will pass this
    ## to the C segmentation code directly.
    counts.by.sample[[sample.id]] <- with(sample.counts, {
      data.frame(chromStart, chromEnd, count)
    })
  }
  counts.by.sample
}

