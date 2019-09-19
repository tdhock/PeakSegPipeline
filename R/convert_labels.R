convert_labels <- function
### Convert label text files in proj.dir/labels/*.txt to
### proj.dir/samples/*/*/labels.bed files.
(proj.dir,
### project directory.
  verbose=0
### Print messages?
){
  chromStart <- chromEnd <- . <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  project.dir <- normalizePath(proj.dir, mustWork=TRUE)
  labels.file.vec <- Sys.glob(file.path(project.dir, "labels", "*.txt"))
  if(length(labels.file.vec)==0){
    stop(
      "proj.dir=", project.dir,
      " should be a project directory with proj.dir/labels/*.txt files")
  }
  label.colors <-
    c(noPeaks="246,244,191",
      peakStart="255,175,175",
      peakEnd="255,76,76",
      peaks="164,69,238")

  sample.dir.glob <- file.path(project.dir, "samples", "*", "*")
  sample.dir.vec <- Sys.glob(sample.dir.glob)
  if(length(sample.dir.vec)==0){
    stop(
      "no sample directories, please create at least one ",
      file.path(
        project.dir, "samples", "groupID", "sampleID", "coverage.bigWig")
    )
  }
  sample.id <- basename(sample.dir.vec)
  sample.group.dir <- dirname(sample.dir.vec)
  sample.group <- basename(sample.group.dir)
  sample.dt <- data.table(sample.id, sample.group)

  regions.by.file <- list()
  chunk.limits.list <- list()
  bed.list <- list()
  positive.regions.list <- list()
  for(labels.file in labels.file.vec){
    if(verbose)cat("Reading ", labels.file, "\n", sep="")
    labels.lines <- readLines(labels.file)
    is.blank <- labels.lines == ""
    chunk.id <- cumsum(is.blank)+1L
    label.df <- data.frame(chunk.id, line=labels.lines)[!is.blank, ]
    if(verbose)cat(length(unique(label.df$chunk.id)), " chunks, ",
        nrow(label.df), " label lines\n", sep="")
    ## Error checking.
    raw.vec <- paste(label.df$line)
    line.vec <- gsub(",", "", raw.vec)
    keep.digits <- function(x)as.integer(gsub("[^0-9]+", "", x))
    int.pattern <- list("[0-9 ,]+", keep.digits)
    match.dt <- nc::capture_first_vec(
      line.vec,
      chrom="chr.+?",
      ":",
      chromStart=int.pattern,
      "-",
      chromEnd=int.pattern,
      " ",
      annotation="[a-zA-Z]+",
      sample_groups=".*")
    not.recognized <- ! match.dt$annotation %in% names(label.colors)
    if(any(not.recognized)){
      print(raw.vec[not.recognized])
      print(match.dt[not.recognized])
      stop(
        "unrecognized annotation; valid values: ",
        paste(
          names(label.colors),
          collapse=", "))
    }
    ## error checking messages should include chunk.id
    match.dt[, chunk.id := label.df$chunk.id ]
    setkey(match.dt, chrom, chromStart, chromEnd)
    match.dt[, {
      same.as.next <- diff(chromStart) <= 0
      if(any(same.as.next)){
        bad.i <- which(same.as.next)
        print(.SD[c(bad.i, bad.i+1), ])
        stop("chromStart not increasing")
      }
      if(any(is.bad <- chromStart >= chromEnd)){
        print(.SD[is.bad])
        stop("chromStart >= chromEnd")
      }
      overlaps.next <- chromStart[-1] < chromEnd[-length(chromEnd)]
      if(any(overlaps.next)){
        print(data.table(.SD, overlaps.next=c(overlaps.next, FALSE)))
        stop("overlapping regions")
      }
    }, by=chrom]
    ## determine total set of sample groups with positive=Peak
    ## annotations.
    stripped <- gsub(" *$", "", gsub("^ *", "", match.dt$sample_groups))
    no.groups <- stripped == ""
    bad.positive <- match.dt$annotation!="noPeaks" & no.groups
    if(any(bad.positive)){
      print(match.dt[bad.positive])
      stop("need at least one sample group up for each positive label")
    }
    commas <- gsub(" +", ",", stripped)
    match.dt[, group.list := strsplit(commas, split=",") ]
    bed.list[[labels.file]] <- match.dt[, data.table(
      chrom,
      chromStart,
      chromEnd,
      name=paste0(annotation, ":", commas),
      score=0,
      strand=".",
      thickStart=chromStart,
      thickEnd=chromEnd,
      itemRgb=label.colors[paste(annotation)])]
    sample.group.vec <- unique(unlist(match.dt$group.list))
    if(verbose)cat("labeled sample groups: ",
        paste(sample.group.vec, collapse=", "),
        "\n",
        sep="")
    labeled.samples <- sample.dt[sample.group.vec, on=.(sample.group)]
    na.samples <- labeled.samples[is.na(sample.id)]
    if(nrow(na.samples)){
      glob.vec <- file.path(na.samples$sample.group, "*")
      glob.str <- paste(glob.vec, collapse=" ")
      stop("no ", glob.str, " directories (but labels are present)")
    }
    ## Create some labeled regions for specific/nonspecific peaks.
    file.positive.regions <- match.dt[0 < sapply(group.list, length)]
    input.has.peak <- grepl("Input", file.positive.regions$sample_groups)
    if(any(input.has.peak)){
      positive.regions.list[[labels.file]] <-
        file.positive.regions[, data.table(
          chrom, regionStart=chromStart, regionEnd=chromEnd,
          annotation, input.has.peak
        )]
    }
    chunkChrom <- paste(match.dt$chrom[1])
    match.dt[, {
      ## Check that all regions are on the same chrom.
      if(any(chrom != chunkChrom)){
        print(.SD)
        stop("each chunk must span only 1 chrom")
      }
    }, by=.(chunk.id)]
    chunk.limits.list[[labels.file]] <- match.dt[, data.table(
      chrom=chunkChrom,
      chromStart=min(chromStart),
      chromEnd=max(chromEnd)
    ), by=.(chunk.id)]
    match.dt[, ann.i := 1:.N]
    regions.by.file[[labels.file]] <- match.dt[, {
      one.chunk <- .SD[, {
        is.observed <- sample.group.vec %in% group.list[[1]]
        group.labels <- rbind(
          if(any(is.observed))data.table(
            sample.group=sample.group.vec[is.observed],
            annotation),
          if(any(!is.observed))data.table(
            sample.group=sample.group.vec[!is.observed],
            annotation="noPeaks"))
        sample.dt[group.labels, data.table(
          sample.id,
          sample.group,
          chrom=chunkChrom,
          chromStart,
          chromEnd,
          annotation),
          on=.(sample.group)]
      }]
    }, by=ann.i]
  }
  all.regions <- do.call(rbind, regions.by.file)

  ## Save positive regions for filtering final peaks.
  positive.regions <- do.call(rbind, positive.regions.list)
  ## positive.regions.RData <- file.path(data.dir, "positive.regions.RData")
  ## save(positive.regions, file=positive.regions.RData)

  ## Save labels to bed file for viewing on UCSC.
  bed <- do.call(rbind, bed.list)[order(chrom, chromStart, chromEnd)]
  all_labels.bed <- file.path(project.dir, "all_labels.bed")
  fwrite(
    bed, all_labels.bed,
    row.names=FALSE,
    col.names=FALSE,
    quote=FALSE)

  ## Write chunk info.
  chunk.limits <- do.call(rbind, chunk.limits.list)
  chunk.limits[, {
    clustered <- clusterPeaks(.SD)
    limits.by.cluster <- split(clustered, clustered$cluster)
    chunks.per.cluster <- sapply(limits.by.cluster, nrow)
    not.ok <- 1 < chunks.per.cluster
    if(any(not.ok)){
      print(limits.by.cluster[not.ok])
      stop("chunks in different label files should not overlap")
    }
  }, by=.(chrom)]
  chunk.limits.csv <- file.path(project.dir, "chunk.limits.csv")
  fwrite(chunk.limits, chunk.limits.csv)

  ## Write labels to each sample.
  all.regions[, {
    labels.bed <- file.path(
      project.dir, "samples", sample.group, sample.id, "labels.bed")
    fwrite(
      data.table(chrom, chromStart, chromEnd, annotation),
      labels.bed,
      quote=FALSE,
      col.names=FALSE,
      sep="\t")
    if(verbose)cat("Wrote ", .N,
        " labels to ", labels.bed,
        "\n", sep="")
  }, by=.(sample.group, sample.id)]

  all.regions
}
