### URL where UCSC data files can be downloaded.
ucsc.goldenPath.url <- "http://hgdownload.soe.ucsc.edu/goldenPath/"

downloadProblems <- function
### Create problems.bed based on data from UCSC.
(genome,
### UCSC ID e.g. hg19 or hg38.
  problems.bed,
### file to save.
  url.prefix=ucsc.goldenPath.url
### http://path.to/goldenPath/
  ){
  gap <- chromInfo <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  pre.genome <- paste0(url.prefix, genome)
  file.list <- list()
  for(db in c("chromInfo", "gap")){
    u <- paste0(pre.genome, "/database/", db, ".txt.gz")
    f <- tempfile()
    gz <- paste0(f, ".gz")
    download.file(u, gz)
    system.or.stop(paste("gunzip", shQuote(gz)))
    if(db=="gap"){
      dt <- fread(f)
      fwrite(dt[, 2:4, with=FALSE], f)
    }
    file.list[[db]] <- f
  }
  with(file.list, gap2problems(gap, chromInfo, problems.bed))
}
