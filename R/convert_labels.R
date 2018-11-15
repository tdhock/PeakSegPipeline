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
  g.pos.pattern <-
    paste0("(?<chrom>chr.+?)",
           ":",
           "(?<chromStart>[0-9 ,]+)",
           "-",
           "(?<chromEnd>[0-9 ,]+)",
           " ",
           "(?<annotation>[a-zA-Z]+)",
           "(?<sample_groups>.*)")
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
  sample.df <- data.frame(sample.id, sample.group)
  samples.by.group <- split(sample.df, sample.df$sample.group)

  regions.by.file <- list()
  regions.by.chunk.file <- list()
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
    match.mat <- namedCapture::str_match_named(line.vec, g.pos.pattern)
    if(any(is.bad.line <- is.na(match.mat[,1]))){
      print(raw.vec[is.bad.line])
      stop("label line does not match ", g.pos.pattern)
    }
    not.recognized <- ! match.mat[, "annotation"] %in% names(label.colors)
    if(any(not.recognized)){
      print(raw.vec[not.recognized])
      print(match.mat[not.recognized, , drop=FALSE])
      stop(
        "unrecognized annotation; valid values: ",
        paste(
          names(label.colors),
          collapse=", "))
    }
    match.df <-
      data.frame(chrom=match.mat[, "chrom"],
                 chromStart=as.integer(match.mat[, "chromStart"]),
                 chromEnd=as.integer(match.mat[, "chromEnd"]),
                 annotation=match.mat[, "annotation"],
                 sample.groups=match.mat[, "sample_groups"],
                 chunk.id=label.df$chunk.id,
                 stringsAsFactors=FALSE)
    match.by.chrom <- split(match.df, match.df$chrom)
    for(chrom in names(match.by.chrom)){
      chrom.df <- match.by.chrom[[chrom]]
      sorted <- chrom.df[with(chrom.df, order(chromStart, chromEnd)), ]
      same.as.next <- diff(sorted$chromStart) <= 0
      if(any(same.as.next)){
        bad.i <- which(same.as.next)
        print(sorted[c(bad.i, bad.i+1), ])
        stop("chromStart not increasing")
      }
      if(any(with(sorted, chromStart >= chromEnd))){
        print(sorted)
        stop("chromStart >= chromEnd")
      }
      overlaps.next <-
        with(sorted, chromStart[-1] < chromEnd[-length(chromEnd)])
      if(any(overlaps.next)){
        print(data.frame(sorted, overlaps.next=c(overlaps.next, FALSE)))
        stop("overlapping regions")
      }
    }

    ## determine total set of sample groups with positive=Peak
    ## annotations.
    stripped <- gsub(" *$", "", gsub("^ *", "", match.df$sample.groups))
    no.groups <- stripped == ""
    bad.positive <- match.df$annotation!="noPeaks" & no.groups
    if(any(bad.positive)){
      print(match.df[bad.positive,])
      stop("need at least one sample group up for each positive label")
    }
    commas <- gsub(" +", ",", stripped)
    sample.group.list <- strsplit(commas, split=",")
    bed.list[[labels.file]] <- 
      data.table(match.df[,c("chrom", "chromStart", "chromEnd")],
                 name=paste0(match.df$annotation, ":", commas),
                 score=0,
                 strand=".",
                 thickStart=match.df$chromStart,
                 thickEnd=match.df$chromEnd,
                 itemRgb=label.colors[paste(match.df$annotation)])
    names(sample.group.list) <- rownames(match.df)
    sample.group.vec <- unique(unlist(sample.group.list))
    if(verbose)cat("labeled sample groups: ",
        paste(sample.group.vec, collapse=", "),
        "\n",
        sep="")
    ## Create some labeled regions for specific/nonspecific peaks.
    groups.up.vec <- sapply(sample.group.list, length)
    file.positive.regions <- match.df[0 < groups.up.vec,]
    input.has.peak <- grepl("Input", file.positive.regions$sample.groups)
    if(any(input.has.peak)){                         
      positive.regions.list[[labels.file]] <-
        with(file.positive.regions, data.frame(
          chrom, regionStart=chromStart, regionEnd=chromEnd,
          annotation, input.has.peak
        ))
    }

    match.by.chunk <- split(match.df, match.df$chunk.id)
    for(chunk.id in names(match.by.chunk)){
      ## Check that all regions are on the same chrom.
      chunk.df <- match.by.chunk[[chunk.id]]
      chunkChrom <- paste(chunk.df$chrom[1])
      if(any(chunk.df$chrom != chunkChrom)){
        print(chunk.df)
        stop("each chunk must span only 1 chrom")
      }
      regions.list <- list()
      for(ann.i in 1:nrow(chunk.df)){
        chunk.row <- chunk.df[ann.i, ]
        groups.up.vec <- sample.group.list[[rownames(chunk.row)]]
        is.observed <- sample.group.vec %in% groups.up.vec
        observed <- sample.group.vec[is.observed]
        not.observed <- sample.group.vec[!is.observed]
        to.assign <- list()
        ann <- chunk.row$annotation
        to.assign[observed] <- ann
        to.assign[not.observed] <- "noPeaks"
        for(sample.group in names(to.assign)){
          relevant.samples <- samples.by.group[[sample.group]]
          if(length(relevant.samples) == 0){
            glob.str <- file.path(sample.group, "*")
            stop("no ", glob.str, " directories (but labels are present)")
          }
          annotation <- to.assign[[sample.group]]
          regions.list[[paste(ann.i, sample.group)]] <- 
            data.table(sample.id=paste(relevant.samples$sample.id),
                       sample.group,
                       chrom=chunk.row$chrom,
                       chromStart=chunk.row$chromStart,
                       chromEnd=chunk.row$chromEnd,
                       annotation)
        }
      }
      one.chunk <- do.call(rbind, regions.list)
      setkey(one.chunk, chromStart, chromEnd)
      file.and.chunk <- paste0(basename(labels.file), "-chunk", chunk.id)
      regions.by.chunk.file[[file.and.chunk]] <- one.chunk
      regions.by.file[[labels.file]][[chunk.id]] <- one.chunk
      chunk.limits.list[[file.and.chunk]] <- with(one.chunk, {
        data.frame(file.and.chunk,
                   chrom=chrom[1],
                   chromStart=min(chromStart),
                   chromEnd=max(chromEnd))
      })
    }
  }
  bed <- do.call(rbind, bed.list)[order(chrom, chromStart, chromEnd),]
  chunk.limits <- do.call(rbind, chunk.limits.list)
  positive.regions <- do.call(rbind, positive.regions.list)
  rownames(chunk.limits) <- NULL
  rownames(bed) <- NULL
  rownames(positive.regions) <- NULL

  ## Save positive regions for filtering final peaks.
  ## positive.regions.RData <- file.path(data.dir, "positive.regions.RData")
  ## save(positive.regions, file=positive.regions.RData)

  ## Save labels to bed file for viewing on UCSC.
  all_labels.bed <- file.path(project.dir, "all_labels.bed")
  ## con <- file(all_labels.bed, "w")
  ## header <- 
  ##   paste("track",
  ##         "visibility=pack",
  ##         "name=PeakSegJointLabels",
  ##         'description="Visually defined labels',
  ##         'in regions with and without peaks"',
  ##         "itemRgb=on")
  ## writeLines(header, con)
  ## write.table(bed, con,
  ##             row.names=FALSE,
  ##             col.names=FALSE,
  ##             quote=FALSE)
  ## close(con)
  write.table(bed, all_labels.bed,
              row.names=FALSE,
              col.names=FALSE,
              quote=FALSE)

  limits.by.chrom <- split(chunk.limits, chunk.limits$chrom)
  for(chrom in names(limits.by.chrom)){
    chrom.limits <- limits.by.chrom[[chrom]]
    ## Find overlapping chunks, and join them:
    clustered <- clusterPeaks(chrom.limits)
    limits.by.cluster <- split(clustered, clustered$cluster)
    chunks.per.cluster <- sapply(limits.by.cluster, nrow)
    not.ok <- 1 < chunks.per.cluster
    if(any(not.ok)){
      print(limits.by.cluster[not.ok])
      stop("chunks in different label files should not overlap")
    }
  }

  ## Write labels to each sample.
  all.regions <- do.call(rbind, regions.by.chunk.file)
  regions.by.sample <- split(all.regions, all.regions[, paste0(sample.group, "/", sample.id)])
  for(sample.path in names(regions.by.sample)){
    sample.labels <- regions.by.sample[[sample.path]]
    labels.bed <- file.path(project.dir, "samples", sample.path, "labels.bed")
    write.table(
      sample.labels[, .(chrom, chromStart, chromEnd, annotation)],
      labels.bed,
      quote=FALSE,
      row.names=FALSE,
      col.names=FALSE,
      sep="\t")
    if(verbose)cat("Wrote ", nrow(sample.labels),
        " labels to ", labels.bed,
        "\n", sep="")
  }

  ## Write chunk info.
  chunk.limits.RData <- file.path(project.dir, "chunk.limits.RData")
  save(chunk.limits, regions.by.chunk.file, file=chunk.limits.RData)
}
