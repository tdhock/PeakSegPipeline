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

  ## Create a trackDb.txt file with links to a few small bigWig files on GitHub.
  bigWig.part.vec <- c(
    "Input/MS010302",
    "bcell/MS010302",
    "Input/MS002201",
    "kidney/MS002201")
  repos.url <- "https://raw.githubusercontent.com/tdhock/input-test-data/master/"
  track.lines <- sprintf(
    "track %s\nbigDataUrl %s%s.bigwig\n",
    sub("/", "_", bigWig.part.vec),
    repos.url, bigWig.part.vec)
  track.dir <- tempfile()
  dir.create(track.dir)
  trackDb.txt <- file.path(track.dir, "trackDb.txt")
  cat(track.lines, sep="\n")
  writeLines(track.lines, trackDb.txt)

  ## Download the bigWig files mentioned in trackDb.txt
  if(interactive()){
    downloadBigWigs(trackDb.txt, track.dir)
  }

})
