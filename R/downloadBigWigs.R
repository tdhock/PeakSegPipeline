downloadBigWigs <- function(trackDb.txt){
  track.pattern <- paste0(
    "    ",
    "(?<name>.*?)",
    " ",
    "(?<value>.*?)",
    "\n")
  trackDb.vec <- readLines(trackDb.txt)
  trackDb.str <- paste(trackDb.vec, collapse="\n")
  track.vec <- strsplit(trackDb.str, "\n\\s*\n")[[1]]
  subtrack.vec <- paste0(track.vec[-1], "\n")
  track.list <- str_match_all_named(subtrack.vec, track.pattern)
  subGroup.pattern <- paste0(
    "(?<name>[^ ]+)",
    "=",
    "(?<value>[^ ]+)")
  data.dir <- dirname(trackDb.txt)
  for(track.i in seq_along(track.list)){
    track.mat <- track.list[[track.i]]
    cat(sprintf("%4s / %4s tracks\n", track.i, length(track.list)))
    if(track.mat["type", "value"] == "bigWig"){
      subGroup.mat <- str_match_all_named(
        track.mat["subGroups", "value"],
        subGroup.pattern)[[1]]
      cell.type <- match["sampleType", "value"]
      u <- track.mat["bigDataUrl", "value"]
      bigWig.base <- sub("[.][^.]+$", ".bigWig", basename(u))
      sample.id <- track.mat["shortLabel", "value"]
      sample.dir <- file.path(data.dir, "samples", cell.type, sample.id)
      dir.create(sample.dir, showWarnings=FALSE)
      bigWig.path <- file.path(sample.dir, "coverage.bigWig")
      if(!file.exists(bigWig.path)){
        dir.create(sample.dir, showWarnings=FALSE, recursive=TRUE)
        download.file(u, bigWig.path)
      }
    }
  }
}
