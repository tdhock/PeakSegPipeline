denormalizeBedGraph <- function
### Attempt to convert a bedGraph file with non-negative non-integer
### data to integer count data, by dividing each value by the smallest
### non-zero value.
(input.bedGraph,
### Input bedGraph file with non-integer data.
  output.bedGraph
### Output bedGraph file with integer data.
){
  stopifnot(is.character(input.bedGraph))
  stopifnot(is.character(output.bedGraph))
  input.path <- normalizePath(input.bedGraph, mustWork=TRUE)
  output.path <- normalizePath(output.bedGraph, mustWork=FALSE)
  result <- .C(
    "denormalize_interface",
    input.path,
    output.path,
    PACKAGE="PeakSegPipeline")
### Nothing.
}

denormalizeBigWig <- function
### Attempt to convert a bigWig file with non-negative non-integer
### data to integer count data, by dividing each value by the smallest
### non-zero value.
(input.bigWig,
### Input bigWig file with non-integer data.
  output.bigWig
### Output bigWig file with integer data.
){
  chrom <- chromEnd <- NULL
  ## Above to avoid CRAN check NOTE.
  if(!(
    is.character(input.bigWig) &&
    length(input.bigWig) == 1 &&
    file.exists(input.bigWig)
  )) {
    stop("input.bigWig must be the name of a data file")
  }
  input.bedGraph <- sub("bigWig$", "bedGraph", input.bigWig)
  stopifnot(is.character(output.bigWig))
  output.bedGraph <- sub("bigWig$", "bedGraph", output.bigWig)
  chromInfo <- bigWigInfo(input.bigWig)
  cmd <- paste(
    "bigWigToBedGraph",
    input.bigWig,
    input.bedGraph)
  system.or.stop(cmd)
  denormalizeBedGraph(input.bedGraph, output.bedGraph)
  cmd <- paste(
    "head",
    input.bedGraph,
    output.bedGraph)
  system.or.stop(cmd)
  output.chromInfo <- sub("bedGraph$", "chromInfo", output.bedGraph)
  chromSizes <- chromInfo[, list(chrom, chromEnd)]
  fwrite(chromSizes, output.chromInfo, sep="\t", col.names=FALSE)
  bedGraphToBigWig(
    output.bedGraph,
    output.chromInfo,
    output.bigWig)
  unlink(input.bedGraph)
  unlink(output.bedGraph)
}
