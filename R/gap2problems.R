gap2problems <- function
(gap.bed,
### gap.bed file, parts of the genome with gaps (NNNN).
  chromInfo.txt,
### chromInfo.txt file, chromosome sizes.
  problems.bed
### will be created
){
  gap.all <- fread(gap.bed)
  gap <- gap.all[, 1:3, with=FALSE]
  setnames(gap, c("chrom", "chromStart", "chromEnd"))
  setkey(gap, chrom)
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
