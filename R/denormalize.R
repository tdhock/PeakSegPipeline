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
  output.bigWig,
### Output bigWig file with integer data.
  verbose=getOption("PeakSegPipeline.verbose", 1)
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
  in.path.bigWig <- normalizePath(input.bigWig, mustWork=TRUE)
  out.path.bigWig <- normalizePath(output.bigWig, mustWork=FALSE)
  input.bedGraph <- sub("bigWig$", "bedGraph", in.path.bigWig)
  stopifnot(is.character(out.path.bigWig))
  output.bedGraph <- sub("bigWig$", "bedGraph", out.path.bigWig)
  chromInfo <- bigWigInfo(in.path.bigWig)
  bigWigToBedGraph(in.path.bigWig, input.bedGraph)
  denormalizeBedGraph(input.bedGraph, output.bedGraph)
  if(verbose)system.or.stop(paste(
    "head",
    shQuote(input.bedGraph),
    shQuote(output.bedGraph)))
  output.chromInfo <- sub("bedGraph$", "chromInfo", output.bedGraph)
  chromSizes <- chromInfo[, list(chrom, chromEnd)]
  fwrite(chromSizes, output.chromInfo, sep="\t", col.names=FALSE)
  bedGraphToBigWig(
    output.bedGraph,
    output.chromInfo,
    out.path.bigWig)
  unlink(input.bedGraph)
  unlink(output.bedGraph)
}
