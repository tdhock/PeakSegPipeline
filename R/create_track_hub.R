create_track_hub <- function
### Create track hub for a project.
(data.dir.path,
### data/project directory path on local filesystem.
  data.dir.url,
### URL for data/project directory, for constructing URLs to
### bigWig/bigBed files in the trackDb.txt file, e.g. if
### data.dir.url="http://some.domain/~user/demo_project" then your
### trackDb.txt file will contain URLs like
### http://some.domain/~user/demo_project/samples/groupID/sampleID/coverage.bigWig
  genome,
### genome string as defined at UCSC, e.g. "hg19"
  email,
### email address for maintainer of track hub.
  goldenPath.url=ucsc.goldenPath.url
### link to download UCSC genome chromInfo files, necessary for
### creating bigWigs.
){
  chrom <- problemStart <- problemEnd <- chromEnd <- . <-
    name <- chrom <- chromStart <- strand <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  chromInfo.dir <- file.path(data.dir.path, "chromInfo")
  dir.create(chromInfo.dir, showWarnings=FALSE)
  chromInfo.txt <- file.path(
    chromInfo.dir,
    paste0(genome, "_chromInfo.txt"))
  ## First make sure we have the chromInfo file for this genome.
  if(!file.exists(chromInfo.txt)){
    chromInfo.url <- paste0(goldenPath.url, genome, "/database/chromInfo.txt.gz")
    chromInfo.gz <- paste0(chromInfo.txt, ".gz")
    download.file(chromInfo.url, chromInfo.gz)
    system.or.stop(paste("zcat", chromInfo.gz, ">", chromInfo.txt))
  }
  ## Then create bedGraph files if necessary.
  bedGraph.file.vec <- Sys.glob(file.path(
    data.dir.path, "samples", "*", "*", "coverage.bedGraph"))
  for(bedGraph.file in bedGraph.file.vec){
    bigWig <- sub("bedGraph$", "bigWig", bedGraph.file)
    if(!file.exists(bigWig)){
      bedGraphToBigWig(bedGraph.file, chromInfo.txt, bigWig)
    }
  }
  bigWig.glob <- file.path(data.dir.path, "samples", "*", "*", "coverage.bigWig")
  bigWig.file.vec <- Sys.glob(bigWig.glob)
  if(length(bigWig.file.vec)==0){
    stop("no ", bigWig.glob, " files")
  }
  getURL <- function(file.vec){
    paste0(data.dir.url, sub(data.dir.path, "", file.vec))
  }
  url.vec <- getURL(bigWig.file.vec)
  sample.path.vec <- dirname(bigWig.file.vec)
  sample.id.vec <- basename(sample.path.vec)
  group.path.vec <- dirname(sample.path.vec)
  group.id.vec <- basename(group.path.vec)
  group.names <- unique(group.id.vec)
  ##dput(RColorBrewer::brewer.pal(Inf, "Set3"))
  maybe.short <- c(
    "#8DD3C7",
    ##"#FFFFB3",#yellow
    "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
    "#B3DE69", "#FCCDE5",
    "#D9D9D9",#grey
    "#BC80BD", "#CCEBC5", "#FFED6F"
  )
  group.colors <- rep(maybe.short, l=length(group.names))
  names(group.colors) <- group.names
  data.name <- basename(data.dir.path)
  joint_peaks.bedGraph.vec <- sub(
    "coverage.bigWig$", "joint_peaks.bedGraph", bigWig.file.vec)
  joint.bigWig.list <- list()
  for(joint_peaks.bedGraph in joint_peaks.bedGraph.vec){
    joint_peaks.bigWig <- sub("bedGraph$", "bigWig", joint_peaks.bedGraph)
    created <- bedGraphToBigWig(
      joint_peaks.bedGraph, chromInfo.txt, joint_peaks.bigWig)
    if(created){
      joint.bigWig.list[[joint_peaks.bedGraph]] <- joint_peaks.bigWig
    }
  }
  ## Write genomes.txt
  writeLines(paste0("
genome ", genome, "
trackDb trackDb.txt
"), file.path(data.dir.path, "genomes.txt"))
  ## Write hub.txt
  writeLines(paste0("
hub ", data.name, "
shortLabel ", data.name, "
longLabel ", data.name, "
genomesFile genomes.txt
email ", email), file.path(data.dir.path, "hub.txt"))
  ## create jointProblems.bigBed
  jproblems.glob <- file.path(data.dir.path, "problems", "*", "jointProblems.bed")
  jprobs <- fread(cmd=paste("cat", jproblems.glob))
  jointProblems.bed <- file.path(data.dir.path, "jointProblems.bed")
  if(nrow(jprobs)){
    setnames(jprobs, c("chrom", "problemStart", "problemEnd"))
    sizes.dt <- fread(chromInfo.txt)
    names(sizes.dt)[1:2] <- c("chrom", "chromEnd")
    join.dt <- sizes.dt[jprobs, on=list(chrom)]
    join.dt[, problemStart := ifelse(problemStart < 0, 0, problemStart)]
    join.dt[, problemEnd := ifelse(problemEnd < chromEnd, problemEnd, chromEnd)]
    setkey(join.dt, chrom, problemStart, problemEnd)
    write.table(
      join.dt[, .(chrom, problemStart, problemEnd)],
      jointProblems.bed,
      quote=FALSE,
      row.names=FALSE,
      col.names=FALSE)
  }else{
    unlink(jointProblems.bed)
  }
  bedToBigBed <- function(bed){
    bed.long <- fread(bed)
    names(bed.long)[1:3] <- c("chrom", "chromStart", "chromEnd")
    if(4 <= ncol(bed.long)){
      names(bed.long)[4] <- "name"
      bed.long[, name := substr(name, 1, 255)]
    }
    if(6 <= ncol(bed.long) && is.numeric(bed.long[[6]])){
      ## a column of all . is read as numeric 0 by fread.
      names(bed.long)[6] <- "strand"
      bed.long[, strand := rep(".", .N)]
    }
    short <- sub(".bed$", "-short.bed", bed)
    setkey(bed.long, chrom, chromStart)
    fwrite(bed.long, short, sep="\t", col.names=FALSE, quote=FALSE)
    bigBed <- sub("bed$", "bigBed", bed)
    cmd <- paste(
      "bedToBigBed",
      short, chromInfo.txt,
      bigBed)
    system.or.stop(cmd)
    bigBed
  }
  bed.num.vec <- c(
    all_labels=9,
    problems=3,
    jointProblems=3,
    peaks_summary=5)
  long.name.vec <- c(
    all_labels="Manually labeled regions with and without peaks",
    problems=paste(
      "Separate problems",
      "(PeakSegFPOP looks for multiple peaks in each region,",
      "independently for each sample)"),
    jointProblems=paste(
      "Joint problems",
      "(PeakSegJoint looks for one common peak in each region,",
      "across all samples)"),
    peaks_summary="Regions with a peak in at least one sample")
  bigBed.list <- list()
  for(bed.name in names(bed.num.vec)){
    bed.file <- file.path(data.dir.path, paste0(bed.name, ".bed"))
    if(file.exists(bed.file)){
      bigBed.list[[bed.name]] <- bedToBigBed(bed.file)
    }
  }
  bed.track.vec <- if(length(bigBed.list)==0){
    ""
  }else{
    paste0("
track ", names(bigBed.list), "
type bigBed ", bed.num.vec[names(bigBed.list)], "
shortLabel _model_", names(bigBed.list), "
longLabel ", long.name.vec[names(bigBed.list)], "
visibility pack
itemRgb ", ifelse(names(bigBed.list)=="all_labels", "on", "off"), "
spectrum ", ifelse(names(bigBed.list)=="peaks_summary", "on", "off"), "
bigDataUrl ", getURL(unlist(bigBed.list)))
  }

  group.track.vec <- paste0("
track ", group.names, "
superTrack on show
shortLabel ", group.names, "
longLabel ", group.names, " ChIP-seq samples
")
  track <- function(url, data.type, color){
    paste0("
  track ", track.id.vec, data.type, "
  bigDataUrl ", url, "
  shortLabel ", track.id.vec, data.type, "
  longLabel ", group.id.vec, " | ", sample.id.vec, " | ", data.type, "
  parent ", track.id.vec, "
  type bigWig
  color ", color, "
")
  }
  track.id.vec <- paste0(group.id.vec, "_", sample.id.vec)
  track.vec <- paste0("
 track ", track.id.vec, "
 parent ", group.id.vec, "
 container multiWig
 type bigWig
 shortLabel ", track.id.vec, "
 longLabel ", group.id.vec, " | ", sample.id.vec, "
 graphType points
 aggregate transparentOverlay
 showSubtrackColorOnUi on
 maxHeightPixels 25:12:8
 visibility full
 autoScale on
", {
  track(
    url.vec,
    "Coverage",
    apply(col2rgb(group.colors[group.id.vec]), 2, paste, collapse=",")
  )
}, {
  if(length(joint.bigWig.list)==0){
    ""
  }else{
    track(
      getURL(joint.bigWig.list),
      "Peaks",
      "0,0,0"
    )
  }
})
  u.group.vec <- unique(group.id.vec)
  equals.vec <- paste0(u.group.vec, "=", u.group.vec)
  track.content <- paste(
    paste(group.track.vec, collapse="\n"),
    paste(bed.track.vec, collapse="\n"),
    paste(track.vec, collapse="\n"),
    sep="\n\n")

  writeLines(track.content, file.path(data.dir.path, "trackDb.txt"))

  cat("Created ", getURL("/hub.txt"), "\n", sep="")
}
