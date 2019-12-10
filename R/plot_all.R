### order chromosomes.
orderChrom <- function(chrom.vec, ...){
  before <- after <- NULL
  ## Above to avoid CRAN check NOTE.
  stopifnot(is.character(chrom.vec))
  value.vec <- unique(chrom.vec)
  chr.dt <- nc::capture_first_vec(
    value.vec,
    "chr",
    before="[^_]+",
    after="_.*", "?")
  ord.vec <- chr.dt[, order(
    suppressWarnings(as.numeric(before)),
    before,
    after)]
  rank.vec <- seq_along(value.vec)
  names(rank.vec) <- value.vec[ord.vec]
  order(rank.vec[chrom.vec], ...)
}

plot_all <- function
### Gather and plot results of peak calling, generate summary web page
### set.dir.arg/index.html. Labeled chunk plots are created in
### parallel via future.apply::future_lapply. If set.dir.arg/hub.sh
### exists it is called at the end of this function in order to
### generate a track hub based on the peak model files -- it should
### contain something like Rscript -e
### 'PeakSegPipeline::create_track_hub(...)'
(set.dir.arg,
### Path/to/data/dir.
  zoom.out.times=10,
## The UCSC links in the HTML tables will be zoomed out from the peak
## this number of times.
  verbose=getOption("PeakSegPipeline.verbose", 1)
){
  jobPeaks <- jprob.name <- sample.loss.diff <- group.loss.diff <-
    Input.up <- zoomPos <- n.groups.up <- str.groups.up <-
      n.groups.down <- str.groups.down <- separate.problem <-
        img <- n.samples.up <- chunk.limits <- chunk.name <-
          chromStart <- chromEnd <- chromStart1 <- NULL
  problemStart <- problem.name <- chrom <- problemEnd <-
    problem.name <- jobPeaks.RData <- peak.name <- peakStart <-
      peakEnd <- means <- peakBases <- samples.prop <- groups <-
        chisq.pvalue <- loss.diff <- sample.path <- mean.str <-
          . <- most.freq.group <- most.freq.prop <- most.next.prop <-
            sample.counts <- least.freq.group <- least.freq.prop <-
              least.next.prop <- least.next.group <- peak <- samples <-
                annotation <- labelStart <- labelEnd <- n.Input <-
                  prop.noPeaks <- FP <- FN <- specificity <- log10.bases <-
                    log10.samples <- n.samples <- log10.samples.i <-
                      log10.bases.i <- log10.bases.start <- log10.bases.end <-
                        log10.samples.start <- log10.samples.end <-
                          log10.bases.label <- log10.bases.mid <-
                            log10.samples.label <- log10.samples.mid <-
                              count <- e <- mid.peakEnd <- min.peakEnd <-
                                max.peakEnd <- chrom.box <- chrom.fac <-
                                  xval <- chunk <- chunk.problem <-
                                    labeled.chunks <- labelStart1 <-
                                      problemStart1 <- labeled.regions <-
                                        problem <- peak.regions <-
                                          peak.samples <- score <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  zoom.factor <- (zoom.out.times-1)/2
  set.dir <- normalizePath(set.dir.arg, mustWork=TRUE)
  ##load all jobPeaks files.
  joint.glob <- file.path(
    set.dir, "jobs", "*")
  jobPeaks.RData.vec <- Sys.glob(file.path(joint.glob, "jobPeaks.RData"))
  if(length(jobPeaks.RData.vec)==0){
    stop(
      "no predicted joint peaks found; to do joint peak prediction run ",
      file.path(joint.glob, "jobPeaks.sh"))
  }
  if(verbose)cat(
    "Reading predicted peaks in",
    length(jobPeaks.RData.vec),
    "jobPeaks.RData files.\n",
    sep=" ")
  ## First pass to figure out background level for each sample.
  summary.dt.list <- list()
  background.list <- list()
  for(job.i in seq_along(jobPeaks.RData.vec)){
    jobPeaks.RData <- jobPeaks.RData.vec[[job.i]]
    load(jobPeaks.RData)
    if(nrow(jobPeaks)){
      jobPeaks[, n.samples := sapply(jobPeaks$background.mean.vec, nrow)]
      not.all.samples <- jobPeaks[n.samples < max(n.samples)]
      if(nrow(not.all.samples)){
        print(not.all.samples)
        stop("some problems do not have all ", max(jobPeaks$n.samples), " samples")
      }
      bkg.mat <- do.call(cbind, jobPeaks$background.mean.vec)
      background.list[[job.i]] <- rowMeans(bkg.mat, na.rm=TRUE)
    }
  }
  mean.background.vec <- rowMeans(do.call(cbind, background.list), na.rm=TRUE)
  out.name.vec <- c(
    sample.loss.diff.vec="likelihood",
    peak.mean.vec="meanCoverage",
    sample.peaks.vec="sample",
    group.peaks.vec="group",
    log10.norm.height="normHeight")
  out.tsv.vec <- structure(
    file.path(set.dir, paste0("peaks_matrix_", out.name.vec, ".tsv.gz")),
    names=names(out.name.vec))
  ## Second pass to save to .tsv.gz matrices.
  conn.list <- list()
  for(out.col in names(out.tsv.vec)){
    out.tsv <- out.tsv.vec[[out.col]]
    conn.list[[out.col]] <- gzfile(out.tsv, "w")
  }
  for(job.i in seq_along(jobPeaks.RData.vec)){
    jobPeaks.RData <- jobPeaks.RData.vec[[job.i]]
    load(jobPeaks.RData)
    if(nrow(jobPeaks)){
      out.mat.list <- list()
      for(out.col in names(out.name.vec)){
        if(out.col %in% names(jobPeaks)){
          out.mat.list[[out.col]] <- do.call(cbind, jobPeaks[[out.col]])
        }
      }
      peak.mean.mat <- do.call(cbind, jobPeaks$peak.mean.vec)
      out.mat.list$log10.norm.height <- log10(peak.mean.mat/mean.background.vec)
      for(out.col in names(conn.list)){
        out.mat <- as.matrix(t(out.mat.list[[out.col]]))
        if(is.logical(out.mat))out.mat[out.mat==TRUE] <- 1L
        out.df <- data.frame(
          peak.name=rownames(out.mat),
          out.mat,
          check.names=FALSE)
        con <- conn.list[[out.col]]
        write.table(
          out.df,
          con,
          quote=FALSE,
          sep=",",
          row.names=FALSE,
          col.names=seek(con)==0)
      }
      isDown <- function(is.up){
        (!is.up) & names(is.up)!="Input"
      }
      job.summary.dt <- jobPeaks[, list(
        problem.name, jprob.name,
        chrom, peakStart, peakEnd, peakBases=peakEnd-peakStart,
        sample.loss.diff, group.loss.diff,
        str.groups.up=apply(out.mat.list$group.peaks.vec, 2, function(is.up){
          paste(names(is.up)[is.up], collapse=",")
        }),
        str.groups.down=apply(out.mat.list$group.peaks.vec, 2, function(is.up){
          paste(names(is.up)[isDown(is.up)], collapse=",")
        }),
        n.groups.down=apply(out.mat.list$group.peaks.vec, 2, function(is.up){
          sum(isDown(is.up))
        }),
        n.groups.up=colSums(out.mat.list$group.peaks.vec),
        n.samples.up=colSums(out.mat.list$sample.peaks.vec)
      )]
      if("Input" %in% rownames(out.mat.list$group.peaks.vec)){
        job.summary.dt[, Input.up := out.mat.list$group.peaks.vec["Input",] ]
      }else{
        job.summary.dt[, Input.up := FALSE ]
      }
      summary.dt.list[[jobPeaks.RData]] <- job.summary.dt
    }#if(nrow(jobPeaks
  }#for(job.i
  sapply(conn.list, close)
  summary.dt <- do.call(rbind, summary.dt.list)
  ## Plot each labeled chunk.
  chunk.limits.csv <- file.path(set.dir, "chunk.limits.csv")
  unsorted.problems <- fread(file=file.path(set.dir, "problems.bed"))
  setnames(unsorted.problems, c("chrom", "problemStart", "problemEnd"))
  unsorted.problems[, problemStart1 := problemStart +1L]
  unsorted.problems[, problem.name := sprintf(
    "%s:%d-%d", chrom, problemStart, problemEnd)]
  setkey(unsorted.problems, chrom, problemStart1, problemEnd)
  chunk.dir.vec <- if(file.exists(chunk.limits.csv)){
    chunks <- fread(file=chunk.limits.csv)
    chunks[, chunk.name := sprintf("%s:%d-%d", chrom, chromStart, chromEnd)]
    chunks[, chromStart1 := chromStart+1L]
    setkey(chunks, chrom, chromStart1, chromEnd)
    chunks.with.problems <- foverlaps(unsorted.problems, chunks, nomatch=0L)
    chunks.with.problems[, file.path(
      set.dir, "problems", problem.name, "chunks", chunk.name)]
  }
  future.apply::future_lapply(chunk.dir.vec, function(chunk.dir){
    PeakSegPipeline::problem.joint.plot(chunk.dir)
  })
  problems <- unsorted.problems[orderChrom(chrom, problemStart),]
  problems[, problem.name := sprintf(
    "%s:%d-%d", chrom, problemStart, problemEnd)]
  problems[, problem.name := factor(problem.name, problem.name)]
  zcat <- function(suffix, ...){
    base.tsv.gz <- paste0(
      "peaks_matrix_", suffix, ".tsv.gz")
    path.tsv.gz <- file.path(set.dir, base.tsv.gz)
    cmd <- paste("zcat", shQuote(path.tsv.gz))
    fread(cmd=cmd, ...)
  }
  ## Read first row and column of samples matrix to get names.
  header.dt <- zcat("sample", nrows=0)
  col.name.vec <- names(header.dt)
  peak.name.dt <- zcat("sample", select=1)
  pos.dt <- problem.table(peak.name.dt$peak.name)
  setnames(pos.dt, c("chrom", "peakStart", "peakEnd"))
  ord.vec <- pos.dt[, orderChrom(chrom, peakStart, peakEnd)]
  pos.ord.dt <- pos.dt[ord.vec]
  for(sample.i in 2:length(col.name.vec)){
    presence.dt <- zcat("sample", select=sample.i)[ord.vec]
    has.peak <- which(presence.dt[[1]]==1)
    mean.vec <- zcat("meanCoverage", select=sample.i)[ord.vec][[1]]
    bg.dt <- data.table(
      pos.ord.dt,
      meanCoverage=sprintf("%.2f", mean.vec)
    )[has.peak]
    out.path <- col.name.vec[[sample.i]]
    out.file <- file.path(
      set.dir, "samples", out.path, "joint_peaks.bedGraph")
    fwrite(
      bg.dt,
      out.file,
      sep="\t",
      col.names=FALSE,
      quote=FALSE)
  }
  specific.html.list <- list()
  group.name.vec <- rownames(jobPeaks$group.peaks.vec[[1]])
  summary.dt[, zoomPos := sprintf(
    "%s:%d-%d", chrom,
    as.integer(peakStart-peakBases*zoom.factor),
    as.integer(peakEnd+peakBases*zoom.factor))]
  for(sg in group.name.vec){
    specific.html.list[[sg]] <- sprintf(
      '<h3>%s</h3>', sg)
    group.peak.list <- list(
      most=summary.dt[
        n.groups.up==1 & str.groups.up==sg][order(-group.loss.diff)],
      least=summary.dt[
        n.groups.down==1 & str.groups.down==sg][order(-group.loss.diff)])
    for(most.or.least in names(group.peak.list)){
      group.peaks <- group.peak.list[[most.or.least]]
      pre.msg <- paste0(
        "<p>",
        nrow(group.peaks),
        " genomic region",
        ifelse(nrow(group.peaks)==1, " was", "s were"))
      msg <- if(most.or.least=="most"){
        paste0(
          pre.msg,
          " predicted to have a common peak in the ",
          sg, " group, and no peaks in other groups.</p>")
      }else{
        paste0(
          pre.msg,
          " predicted to have no peaks in the ",
          sg, " group, with a peak in all other ",
          "groups.</p>")
      }
      tab <- if(nrow(group.peaks)){
        group.peaks[, peak := sprintf('
<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?position=%s">%s</a>
', zoomPos, zoomPos)]
        xt <- xtable(group.peaks[, .(peak, peakBases, group.loss.diff)])
        print(
          xt, type="html",
          print.results=FALSE,
          sanitize.text.function=identity)
      }else ""
      specific.html.list[[paste(sg, most.or.least)]] <- c(msg, tab)
    }
  }
  specific.html.vec <- do.call(c, specific.html.list)
  figure.png.vec <- Sys.glob(file.path(
    set.dir, "problems", "*", "chunks", "*", "figure-predictions.png"))
  if(0 == length(figure.png.vec)){
    chunk.info <- data.table()
    chunks.html <- ""
  }else{
    relative.vec <- sub("/", "", sub(set.dir, "", figure.png.vec))
    chunk.dir.vec <- dirname(figure.png.vec)
    chunks.dir.vec <- dirname(chunk.dir.vec)
    separate.prob.dir.vec <- dirname(chunks.dir.vec)
    g.pos.pattern <- paste0(
      "(?<chrom>chr.+?)",
      ":",
      "(?<chromStart>[0-9 ,]+)",
      "-",
      "(?<chromEnd>[0-9 ,]+)")
    pos2df <- function(path.vec){
      problem <- basename(path.vec)
      int.pattern <- list("[0-9]+", as.integer)
      dt <- nc::capture_first_vec(
        problem,
        chrom="chr.+?",
        ":",
        chromStart=int.pattern,
        "-",
        chromEnd=int.pattern)
      data.frame(problem, dt)
    }
    chunk.info <- data.table(
      separate=pos2df(separate.prob.dir.vec),
      chunk=pos2df(chunk.dir.vec))
    chunk.info[, problem.name := separate.problem]
    chunk.info[, img := sprintf('
<a href="%s">
  <img src="%s" />
</a>
', relative.vec, sub(".png$", "-thumb.png", relative.vec))]
    chunk.info[, chunk := sprintf({
      '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?position=%s">%s</a>'
    }, chunk.problem, chunk.problem)]
    chunk.counts <- chunk.info[, list(chunks=.N), by=problem.name]
    problems[, labeled.chunks := 0L]
    setkey(chunk.counts, problem.name)
    setkey(problems, problem.name)
    problems[chunk.counts, labeled.chunks := chunk.counts$chunks]
    chunks.xt <- xtable(chunk.info[, .(chunk, img)])
    chunks.html <- print(chunks.xt, type="html", print.results=FALSE, sanitize.text.function=identity)
  }
  labels.bed.vec <- Sys.glob(file.path(
    set.dir, "samples", "*", "*", "labels.bed"))
  all.labels.list <- list()
  for(labels.bed in labels.bed.vec){
    sample.dir <- dirname(labels.bed)
    sample.id <- basename(sample.dir)
    group.dir <- dirname(sample.dir)
    sample.group <- basename(group.dir)
    sample.labels <- fread(file=labels.bed)
    setnames(sample.labels, c("chrom", "labelStart", "labelEnd", "annotation"))
    all.labels.list[[labels.bed]] <- data.table(
      sample.id, sample.group, sample.labels)
  }
  all.labels <- do.call(rbind, all.labels.list)
  label.counts <- all.labels[, list(
    samples=.N
  ), by=.(chrom, labelStart, labelEnd)]
  label.counts[, labelStart1 := labelStart + 1L]
  problems[, problemStart1 := problemStart + 1L]
  setkey(label.counts, chrom, labelStart1, labelEnd)
  setkey(problems, chrom, problemStart1, problemEnd)
  labels.in.problems <- foverlaps(label.counts, problems, nomatch=0L)
  label.problem.counts <- labels.in.problems[, list(
    labels=sum(samples),
    labeled.regions=.N
  ), by=problem.name]
  problems[, labels := 0L]
  setkey(label.problem.counts, problem.name)
  setkey(problems, problem.name)
  problems[label.problem.counts, labels := label.problem.counts$labels]
  problems[, labeled.regions := 0L]
  problems[label.problem.counts, labeled.regions := label.problem.counts$labeled.regions]
  problems[, problem := sprintf({
    '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?position=%s">%s</a>'
  }, problem.name, problem.name)]
  ## add peaks, peak.samples, samples
  separate.counts <- summary.dt[, list(
    peak.regions=.N,
    peak.samples=sum(n.samples.up)
  ), by=list(problem.name)]
  setkey(separate.counts, problem.name)
  problems.out <- separate.counts[problems]
  problems.xt <- xtable(problems.out[, list(
    problem,
    peak.regions, peak.samples,
    labeled.chunks, labeled.regions, labels)])
  problems.html <- print(
    problems.xt, type="html",
    print.results=FALSE,
    sanitize.text.function=identity)
  labels.dir <- file.path(set.dir, "labels")
  labels.txt.vec <- Sys.glob(file.path(labels.dir, "*.txt"))
  labels.rel.vec <- sub(labels.dir, "labels", labels.txt.vec)
  labels.href.vec <- sprintf(
    '<a href="%s">%s</a>', labels.rel.vec, labels.rel.vec)
  html.vec <- c(
    '<title>PeakSegFPOP + PeakSegJoint predictions</title>',
    '<h2>Output: predicted peaks</h2>',
    '<ul>',
    sprintf('<li>
%d genomic regions were found to have
at least one sample with a peak
(<a href="hub.txt">track hub</a>,
 <a href="peaks_summary.tsv">peaks_summary.tsv</a>).
</li>', nrow(summary.dt)),
sprintf('
<li>%d peaks were detected overall, across %d samples.</li>
<li>Peak presence/absence matrices:
 <a href="peaks_matrix_sample.tsv.gz">(regions x samples)</a>,
 <a href="peaks_matrix_group.tsv.gz">(regions x groups)</a>,
</li>
<li><a href="peaks_matrix_meanCoverage.tsv.gz">Peak raw height matrix (regions x samples)</a>, mean coverage in peak, which is not comparable between samples, but is useful for plotting the segment mean on top of the raw coverage/count data,</li>
<li><a href="peaks_matrix_normHeight.tsv.gz">Peak normalized height matrix (regions x samples)</a>, log10(mean coverage in peak/mean coverage in background on the same sample), which is comparable between samples,</li>
<li><a href="peaks_matrix_likelihood.tsv.gz">Peak log-likelihood increase matrix (regions x samples)</a> -- useful for ranking peaks (e.g. taking top 100 in each sample) -- wider and taller peaks are more likely, and have larger values in this matrix -- technically the values are computed as: nll_1 - nll_3, where nll_k is the negative log likelihood of the model with k segments (nll_1 is the model with one segment, no peak; nll_3 is the model with a peak, three segments, one change up and one change down).</li>
',
sum(summary.dt$n.samples.up),
length(col.name.vec)-1),
'</ul>',
'<h2>Input: labeled genomic windows</h2>',
sprintf(
  '%d labeled genomic windows (chunks) were defined in %d label file(s):',
  nrow(chunk.info),
  length(labels.txt.vec)),
labels.href.vec,
chunks.html,
'<h2>Input: genomic segmentation problems</h2>',
sprintf(
  '%d problems were defined in <a href="problems.bed">problems.bed</a>',
  nrow(problems)),
problems.html,
'<h2>Output: predicted group-specific peaks</h2>',
'<p>These are great candidates for re-labeling.</p>',
specific.html.vec
)
  writeLines(html.vec, file.path(set.dir, "index.html"))
  ## Write peaks_summary.bed
  fwrite(
    summary.dt,
    file.path(set.dir, "peaks_summary.tsv"),
    sep="\t")
  ## summary
  peaks.bed <- file.path(set.dir, "peaks_summary.bed")
  bed.dt <- summary.dt[Input.up==FALSE]
  max.samples <- max(bed.dt$n.samples.up)
  bed.dt[, score := as.integer((n.samples.up/max.samples)*1000) ]
  bed.dt[str.groups.up=="", str.groups.up := paste0(
    n.samples.up, "sample", ifelse(n.samples.up==1, "", "s"))]
  fwrite(
    bed.dt[, .(
      chrom, peakStart, peakEnd,
      substr(str.groups.up, 1, 255),#max characters in bigBed name.
      score)],
    peaks.bed,
    sep="\t",
    col.names=FALSE)
  ## finally, create track hub if hub.sh exists.
  hub.sh <- file.path(set.dir, "hub.sh")
  if(file.exists(hub.sh)){
    system.or.stop(paste("bash", shQuote(hub.sh)))
  }
### Nothing.
}
