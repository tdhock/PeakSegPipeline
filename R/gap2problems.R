gap2problems <- function
### Convert a gap file from UCSC to a problems.bed file, which is used
### in the PeakSegPipeline to determine where PeakSegFPOP should be
### run to do peak calling for each sample independently.
(gap.bed,
### gap.bed file, parts of the genome with gaps (NNNN).
  chromInfo.txt,
### chromInfo.txt file, chromosome sizes.
  problems.bed
### will be created
){
  chrom <- chromStart <- chromEnd <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  stopifnot(is.character(gap.bed), length(gap.bed)==1, file.exists(gap.bed))
  gap.all <- fread(gap.bed)
  gap <- gap.all[, 1:3, with=FALSE]
  setnames(gap, c("chrom", "chromStart", "chromEnd"))
  setkey(gap, chrom)
  stopifnot(
    is.character(chromInfo.txt),
    length(chromInfo.txt)==1,
    file.exists(chromInfo.txt))
  chromInfo <- fread(chromInfo.txt)
  chromSizes <- chromInfo[, 1:2, with=FALSE]
  setnames(chromSizes, c("chrom", "bases"))
  setkey(chromSizes, chrom)
  problems <- gap[, {
    bases <- chromSizes[chrom]$bases
    problemStart <- c(0, chromEnd)
    problemEnd <- c(chromStart, bases)
    data.table(problemStart, problemEnd)[problemStart < problemEnd,]
  }, by=chrom]
  out <- problems[, list(
    chrom,
    sprintf("%d", problemStart),
    sprintf("%d", problemEnd)
  )]
  write.table(out, problems.bed, sep="\t", quote=FALSE,
              row.names=FALSE, col.names=FALSE)
}
