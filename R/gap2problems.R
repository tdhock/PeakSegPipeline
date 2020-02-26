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
  chrom <- chromStart <- chromEnd <- gapStart <- gapEnd <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  stopifnot(is.character(gap.bed), length(gap.bed)==1, file.exists(gap.bed))
  gap <- fread(file=gap.bed, select=1:3, col.names=c("chrom", "gapStart", "gapEnd"))
  stopifnot(
    is.character(chromInfo.txt),
    length(chromInfo.txt)==1,
    file.exists(chromInfo.txt))
  chromInfo <- fread(file=chromInfo.txt, select=1:2, col.names=c("chrom", "bases"))
  join.dt <- gap[chromInfo, on=list(chrom)]
  problems <- rbind(
    join.dt[is.na(gapStart), {
      data.table(chrom, problemStart=0, problemEnd=bases)
    }],
    join.dt[!is.na(gapStart), {
      bases <- bases[1]
      problemStart <- c(0, gapEnd)
      problemEnd <- c(gapStart, bases)
      data.table(problemStart, problemEnd)[problemStart < problemEnd,]
    }, by=list(chrom)])
  out <- problems[, list(
    chrom,
    sprintf("%d", problemStart),
    sprintf("%d", problemEnd)
  )]
  write.table(out, problems.bed, sep="\t", quote=FALSE,
              row.names=FALSE, col.names=FALSE)
}
