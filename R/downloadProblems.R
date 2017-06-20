arg.vec <- c("hg19", "hg19_download.bed")
arg.vec <- c("hg38", "hg38_download.bed")
arg.vec <- commandArgs(trailingOnly=TRUE)
if(length(arg.vec) != 2){
  stop("usage: Rscript downloadProblems.R ucsc_genome_id problems.bed")
}
pre <- "http://hgdownload.soe.ucsc.edu/goldenPath/"
genome <- arg.vec[1]
pre.genome <- paste0(pre, genome)

library(data.table)

system.or.stop <- function(cmd){
  cat(cmd, "\n")
  code <- system(cmd)
  if(code != 0){
    stop("non-zero exit code ", code)
  }
}

options(warn=2)

file.list <- list()
for(db in c("chromInfo", "gap")){    
  u <- paste0(pre.genome, "/database/", db, ".txt.gz")
  f <- tempfile()
  gz <- paste0(f, ".gz")
  download.file(u, gz)
  system.or.stop(paste("gunzip", gz))
  if(db=="gap"){
    dt <- fread(f)
    fwrite(dt[, 2:4, with=FALSE], f)
  }
  file.list[[db]] <- f
}

cmd <- with(file.list, paste(
  "Rscript gap2problems.R",
  gap, chromInfo, arg.vec[2]))
system.or.stop(cmd)

