downloadBigWigs <- structure(function
### Download bigWig files from a trackDb file.
(trackDb.txt,
### trackDb text file.
  out.dir
### Output directory.
){
  trackDb.vec <- readLines(trackDb.txt)
  track.mat <- namedCapture::str_match_all_variable(
    trackDb.vec,
    "track ",
    name="[^\n]+",
    "(?:\n[^\n]+)*",
    "\\s+bigDataUrl ",
    bigDataUrl="[^\n]+")
  for(track.i in seq_along(track.mat[, "bigDataUrl"])){
    track.name <- rownames(track.mat)[[track.i]]
    u <- track.mat[track.i, "bigDataUrl"]
    dest <- file.path(out.dir, track.name, basename(u))
    cat(sprintf(
      "%4s / %4s tracks %s -> %s\n",
      track.i, nrow(track.mat),
      u, dest))
    dir.create(dirname(dest), showWarnings=FALSE, recursive=TRUE)
    tryCatch({
      download.file(u, dest)
    }, error=function(e){
      message(e)
    })
  }
### Nothing.
}, ex=function(){

  trackDb.txt <- "https://rcdata.nau.edu/genomic-ml/PeakSegFPOP/labels/H3K36me3_TDH_immune/trackDb.txt"
  downloadBigWigs(trackDb.txt, tempfile())

})
