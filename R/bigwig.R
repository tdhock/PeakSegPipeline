bedGraphToBigWig <- function
### Attempt to create bigWig file.
(bedGraph,
### bedGraph file (input), four tab-separated columns: chrom,
### chromStart, chromEnd, numeric data.
  chromInfo,
### chromInfo file (input), two tab-separated columns: chrom, size in
### bases.
  bigWig
### bigWig file (output).
){
  unlink(bigWig)
  bedGraph.dt <- fread(
    bedGraph,
    col.names=c("chrom", "chromStart", "chromEnd", "value"))
  ## This if statement is needed because bedGraphToBigWig stops with
  ## an error code if there are no data, but we don't want to stop,
  ## even if there are no predicted peaks in some samples.
  if(nrow(bedGraph.dt)){
    sorted.dt <- bedGraph.dt[order(chrom, chromStart)]
    bedGraph.sorted <- paste0(basename(bedGraph), ".sorted")
    fwrite(sorted.dt, bedGraph.sorted, sep="\t", col.names=FALSE, quote=FALSE)
    cmd <- paste(
      "bedGraphToBigWig", bedGraph.sorted,
      chromInfo, bigWig)
    system.or.stop(cmd)
  }
  file.exists(bigWig)
### TRUE if bigWig was created.
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
  cmd <- sprintf(
    "bigWigToBedGraph -chrom=%s -start=%d -end=%d %s /dev/stdout",
    chrom, start, end,
    bigwig.file)
  bg <- fread(cmd=cmd, drop=1)
  if(nrow(bg)==0){
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
  chroms <- fread(cmd=cmd, header=FALSE, sep=" ")
  setnames(chroms, c("chrom", "chrom.int", "chromEnd"))
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

stop.without.ucsc <- function
### Stop with an error if UCSC command line programs are not
### available.
(prog.vec=c(
  "bigWigInfo", "bigWigToBedGraph",
  "bedToBigBed", "bedGraphToBigWig")
### Vector of command line programs to test via base::Sys.which.
){
  path.vec <- Sys.which(prog.vec)
  if(any(not.found <- path.vec=="")){
    stop(
      "UCSC command line programs not found, or not executable: ",
      paste(prog.vec[not.found], collapse=", "))
  }
  path.vec
### Full paths to programs.
}
