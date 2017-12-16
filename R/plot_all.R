plot_all <- function
### Plot results of peak calling, generate summary web page.
(set.dir.arg,
### Path/to/data/dir.
  zoom.out.times=10
## The UCSC links in the HTML tables will be zoomed out from the peak
## this number of times.
){
  problemStart <- problem.name <- chrom <- problemEnd <-
    separate.problem <- jobPeaks.RData <- peak.name <- peakStart <-
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
  orderChrom <- function(chrom.vec, ...){
    stopifnot(is.character(chrom.vec))
    chr.pattern <- paste0(
      "chr",
      "(?<before>[^_]+)",
      "(?<after>_.*)?")
    value.vec <- unique(chrom.vec)
    chr.mat <- namedCapture::str_match_named(value.vec, chr.pattern)
    did.not.match <- is.na(chr.mat[, 1])
    if(any(did.not.match)){
      print(value.vec[did.not.match])
      stop("chroms did not match ", chr.pattern)
    }
    ord.vec <- order(
      suppressWarnings(as.numeric(chr.mat[, "before"])),
      chr.mat[, "before"],
      chr.mat[, "after"])
    rank.vec <- seq_along(value.vec)
    names(rank.vec) <- value.vec[ord.vec]
    order(rank.vec[chrom.vec], ...)
  } 
  factorChrom <- function(chrom.vec){
    u.vec <- unique(chrom.vec)
    ord.vec <- u.vec[orderChrom(u.vec)]
    factor(chrom.vec, ord.vec)
  }
  ##load all jobPeaks files.
  joint.glob <- file.path(
    set.dir, "jobs", "*")
  jobPeaks.RData.dt <- data.table(
    jobPeaks.RData=Sys.glob(file.path(joint.glob, "jobPeaks.RData")))
  if(nrow(jobPeaks.RData.dt)==0){
    stop(
      "no predicted joint peaks found; to do joint peak prediction run ",
      file.path(joint.glob, "jobPeaks.sh"))
  }
  cat(
    "Reading predicted peaks in",
    nrow(jobPeaks.RData.dt),
    "jobPeaks.RData files.\n",
    sep=" ")
  jobPeaks <- jobPeaks.RData.dt[, {
    load(jobPeaks.RData)
    jobPeaks
  }, by=jobPeaks.RData]
  jobPeaks[, peak.name := sprintf("%s:%d-%d", chrom, peakStart, peakEnd)]
  ## height
  nrow.vec <- sapply(jobPeaks$background.peak.means, nrow)
  longest.i <- which(nrow.vec==max(nrow.vec))[1]
  one.mat <- jobPeaks$background.peak.means[[longest.i]]
  na.vec <- rep(NA, nrow(one.mat))
  names(na.vec) <- rownames(one.mat)
  getOrNA <- function(m, col.name){
    na.vec[rownames(m)] <- m[, col.name]
    na.vec
  }
  background.mat <- sapply(jobPeaks$background.peak.means, getOrNA, "background")
  peak.mat <- sapply(jobPeaks$background.peak.means, getOrNA, "peak")
  mean.background.vec <- rowMeans(background.mat, na.rm=TRUE)
  log10.peak.height.mat <- log10(peak.mat / mean.background.vec)
  log10.peak.height.mat[is.na(log10.peak.height.mat)] <- 0
  height.dt <- data.table(
    peak.name=jobPeaks$peak.name,
    t(log10.peak.height.mat))
  peaks_matrix_height.tsv <- file.path(set.dir, "peaks_matrix_height.tsv")
  cat(
    "Writing peak height matrix to",
    peaks_matrix_height.tsv,
    "(",
    nrow(height.dt),
    "peaks x",
    ncol(height.dt)-1,
    "samples",
    ")\n",
    sep=" ")
  fwrite(
    height.dt,
    peaks_matrix_height.tsv,
    sep="\t")
  ## Plot each labeled chunk.
  chunk.dir.vec <- Sys.glob(file.path(
    set.dir, "problems", "*", "chunks", "*"))
  LAPPLY <- lapply
  LAPPLY <- mclapply.or.stop
  LAPPLY(chunk.dir.vec, function(chunk.dir){
    PeakSegPipeline::problem.joint.plot(chunk.dir)
  })
  unsorted.problems <- fread(file.path(set.dir, "problems.bed"))
  setnames(unsorted.problems, c("chrom", "problemStart", "problemEnd"))
  chr.pattern <- paste0(
    "chr",
    "(?<before>[^_]+)",
    "(?<after>_.*)?")
  chr.mat <- namedCapture::str_match_named(unsorted.problems$chrom, chr.pattern)
  problems <- unsorted.problems[order(
    suppressWarnings(as.numeric(chr.mat[, "before"])),
    chr.mat[, "before"],
    chr.mat[, "after"],
    problemStart),]
  problems[, problem.name := sprintf(
    "%s:%d-%d", chrom, problemStart, problemEnd)]
  problems[, separate.problem := factor(problem.name, problem.name)]
  ## Count total samples using directories.
  sample.dir.vec <- Sys.glob(file.path(set.dir, "samples", "*", "*"))
  sample.path.vec <- sub(".*samples/", "", sample.dir.vec)
  sample.group.tab <- table(sub("/.*", "", sample.path.vec))
  sample.group.totals <- data.table(
    sample.group=names(sample.group.tab),
    samples.total=as.integer(sample.group.tab))
  setkey(sample.group.totals, sample.group)
  joint.peaks.dt <- jobPeaks[, {
    mean.vec <- means[[1]]
    list(
      separate.problem=problem.name,
      sample.path=names(mean.vec),
      mean=as.double(mean.vec),
      sample.id=sub(".*/", "", names(mean.vec)),
      sample.group=sub("/.*", "", names(mean.vec)),
      peakBases=peakEnd-peakStart)
  }, by=list(chrom, peakStart, peakEnd, peak.name)]
  group.counts <- joint.peaks.dt[, {
    tab <- sort(table(sample.group))
    nonzero <- tab[tab!=0]
    list(
      n.Input=sum(sample.group=="Input"),
      n.samples=.N,
      n.groups=length(tab),
      sample.counts=paste(paste0(
        names(nonzero), ":", nonzero), collapse=","))
  }, by=list(chrom, peakStart, peakEnd, peak.name, peakBases)]
  group.counts.wide <- dcast(
    joint.peaks.dt, chrom + peakStart + peakEnd + peak.name ~ sample.group, length,
    value.var="peakBases")#to avoid message.
  group.counts.mat <- as.matrix(
    group.counts.wide[, sample.group.totals$sample.group, with=FALSE])
  rownames(group.counts.mat) <- group.counts.wide$peak.name
  group.prop.mat <- group.counts.mat / matrix(
    sample.group.totals$samples.total,
    nrow(group.counts.mat),
    ncol(group.counts.mat),
    byrow=TRUE)
  group.prop.tall <- melt(data.table(
    group.counts.wide[,list(chrom, peakStart, peakEnd)],
    group.prop.mat),
    id.vars=c("chrom", "peakStart", "peakEnd"),
    variable.name="sample.group",
    value.name="samples.prop")
  setkey(group.prop.tall, chrom, peakStart, peakEnd, samples.prop)
  group.prop.groups <- group.prop.tall[, list(
    groups=paste(sample.group, collapse=",")
  ), by=list(chrom, peakStart, peakEnd, samples.prop)]
  most.least <- group.prop.groups[, data.table(
    most.freq.group=groups[.N],
    most.freq.prop=samples.prop[.N],
    most.next.group=if(.N==1)NA_character_ else groups[.N-1],
    most.next.prop=if(.N==1)NA_real_ else samples.prop[.N-1],
    least.freq.group=groups[1],
    least.freq.prop=samples.prop[1],
    least.next.group=if(.N==1)NA_character_ else groups[2],
    least.next.prop=if(.N==1)NA_real_ else samples.prop[2]
  ), by=list(chrom, peakStart, peakEnd)]
  input.pred <- most.least[group.counts, on=list(chrom, peakStart, peakEnd)]
  getPvalue <- function(samples.with.peaks, fun){
    count.mat <- rbind(
      samples.with.peaks,
      sample.group.totals$samples.total-samples.with.peaks)
    tryCatch({
      exact <- suppressWarnings(fun(count.mat))
      exact$p.value
    }, error=function(e){
      NA_real_
    })
  }
  if(FALSE){
    input.pred[, fisher.pvalue := apply(
      group.counts.mat[peak.name,], 1, getPvalue, fisher.test)]
  }
  input.pred[, chisq.pvalue := apply(
    group.counts.mat[peak.name,], 1, getPvalue, chisq.test)]
  setkey(jobPeaks, peak.name)
  input.pred[, loss.diff := jobPeaks[input.pred$peak.name, loss.diff] ]
  input.pred[, separate.problem := jobPeaks[input.pred$peak.name, problem.name] ]
  input.pred[, peakBases := peakEnd - peakStart]
  peak.mat <- matrix(
    FALSE,
    length(sample.path.vec),
    nrow(input.pred),
    dimnames=list(
      sample=sample.path.vec,
      peak=input.pred$peak.name))
  i.mat <- joint.peaks.dt[, cbind(sample.path, peak.name)]
  peak.mat[i.mat] <- TRUE
  ## Save samples/groupID/sampleID/joint_peaks.bedGraph files.
  setkey(joint.peaks.dt, sample.path, chrom, peakStart, peakEnd)
  out.path.vec <- unique(joint.peaks.dt$sample.path)
  joint.peaks.dt[, mean.str := sprintf("%.2f", mean)]
  for(out.path in out.path.vec){
    sample.peaks <- joint.peaks.dt[out.path]
    out.file <- file.path(set.dir, "samples", out.path, "joint_peaks.bedGraph")
    fwrite(
      sample.peaks[, .(chrom, peakStart, peakEnd, mean.str)],
      out.file,
      sep="\t",
      col.names=FALSE,
      quote=FALSE)
  }
  specific.html.list <- list()
  for(sg in sample.group.totals$sample.group){
    specific.html.list[[sg]] <- sprintf(
      '<h3>%s</h3>', sg)
    group.peak.list <- list(
      most=input.pred[
        most.freq.group==sg &
        most.freq.prop==1 &
        most.next.prop==0,][order(-loss.diff), .(
          peak.name, loss.diff, samples=sample.counts,
          chrom, peakStart, peakEnd, peakBases)],
      least=input.pred[
        least.freq.group==sg &
        least.freq.prop==0 &
        least.next.prop==1,][order(-loss.diff), .(
          peak.name, loss.diff, samples=least.next.group,
          chrom, peakStart, peakEnd, peakBases)])
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
          " predicted to have a peak in each ",
          sg, " sample, and no peaks in other samples.</p>")
      }else{
        paste0(
          pre.msg,
          " predicted to have no peaks in any ",
          sg, " samples, and at least one other",
          " group with peaks in all samples.</p>")
      }
      tab <- if(nrow(group.peaks)){
        group.peaks[, peak := sprintf('
<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?position=%s:%d-%d">%s</a>
', chrom,
as.integer(peakStart-peakBases*zoom.factor),
as.integer(peakEnd+peakBases*zoom.factor),
peak.name)]
        xt <- xtable(group.peaks[, .(peak, peakBases, loss.diff, samples)])
        print(
          xt, type="html", sanitize.text.function=identity)
      }else ""
      specific.html.list[[paste(sg, most.or.least)]] <- c(msg, tab)   
    }
  }
  specific.html.vec <- do.call(c, specific.html.list)
  labels.bed.vec <- Sys.glob(file.path(
    set.dir, "samples", "*", "*", "labels.bed"))
  all.labels.list <- list()
  for(labels.bed in labels.bed.vec){
    sample.dir <- dirname(labels.bed)
    sample.id <- basename(sample.dir)
    group.dir <- dirname(sample.dir)
    sample.group <- basename(group.dir)
    sample.labels <- fread(labels.bed)
    setnames(sample.labels, c("chrom", "labelStart", "labelEnd", "annotation"))
    all.labels.list[[labels.bed]] <- data.table(
      sample.id, sample.group, sample.labels)
  }
  all.labels <- do.call(rbind, all.labels.list)
  input.labels <- all.labels[sample.group=="Input", list(
    prop.noPeaks=mean(annotation=="noPeaks")
  ), by=.(chrom, labelStart, labelEnd)]
  setkey(input.labels, chrom, labelStart, labelEnd)
  if(nrow(input.labels)){
    setkey(input.pred, chrom, peakStart, peakEnd)
    labeled.input <- foverlaps(input.pred, input.labels, nomatch=0L)
    tryCatch({
      thresh.dt <- labeled.input[, data.table(WeightedROC::WeightedROC(
        n.Input, ifelse(prop.noPeaks==0, 1, -1)))]
      thresh.best <- thresh.dt[which.min(FP+FN),]
      ## threshold is smallest n.Input that is classified as non-specific.
      input.pred[, specificity := ifelse(
        n.Input >= thresh.best$threshold, "non-specific", "specific")]
    }, error=function(e){
      input.pred[, specificity := "unknown"]
    })
  }else{
    input.pred[, specificity := "unknown"]
  }
  edge <- function(x){
    r <- range(input.pred[[x]])
    l <- 36
    e <- seq(r[1], r[2], l=l)
    L <- list()
    L[[paste0(x, ".i")]] <- 1:(l-1)
    start <- e[-l]
    L[[paste0(x, ".start")]] <- start
    end <- e[-1]
    L[[paste0(x, ".end")]] <- end
    mid <- (start+end)/2
    L[[paste0(x, ".mid")]] <- mid
    L[[paste0(x, ".label")]] <- 10^mid
    ##sprintf("%d-%d", ceiling(10^start), floor(10^end))
    do.call(data.table, L)
  }
  input.pred[, log10.bases := log10(peakBases)]
  input.pred[, log10.samples := log10(n.samples)]
  bases.dt <- edge("log10.bases")
  samples.dt <- edge("log10.samples")
  samples.dt
  grid.dt <- data.table(expand.grid(
    log10.bases.i=1:nrow(bases.dt),
    log10.samples.i=1:nrow(samples.dt)))[samples.dt, on=list(
      log10.samples.i)][bases.dt, on=list(log10.bases.i)]
  join.dt <- grid.dt[input.pred, on=list(
    log10.bases.start <= log10.bases,
    log10.bases.end >= log10.bases,
    log10.samples.start <= log10.samples,
    log10.samples.end >= log10.samples)]
  stopifnot(nrow(join.dt)==nrow(input.pred))
  join.counts <- join.dt[, list(
    count=.N
  ), by=list(
    log10.bases.label, log10.bases.mid,
    log10.samples.label, log10.samples.mid)]
  join.biggest <- join.counts[, {
    .SD[which.max(count)]
  }, by=list(log10.samples.mid)]
  chrom.boxes <- input.pred[, {
    m <- max(peakEnd)
    box.size <- 1e7
    n.boxes <- ceiling(m/box.size)
    data.table(e=(0:n.boxes)*box.size)[, data.table(
      min.peakEnd=e[-.N],
      max.peakEnd=e[-1]
    )]
  }, by=list(chrom)]
  chrom.boxes[, mid.peakEnd := (min.peakEnd+max.peakEnd)/2]
  chrom.boxes[, chrom.box := paste0(
    chrom, ":", scales::comma(min.peakEnd), "-", scales::comma(max.peakEnd))]
  peak.boxes <- chrom.boxes[input.pred, on=list(
    chrom,
    min.peakEnd <= peakEnd,
    max.peakEnd >= peakEnd)]
  stopifnot(nrow(peak.boxes)==nrow(input.pred))
  peak.box.counts <- peak.boxes[, list(
    count=.N,
    some.input=sum(0 < n.Input)
  ), by=list(chrom, mid.peakEnd, chrom.box)]
  peak.box.counts[, chrom.fac := factorChrom(chrom)]
  getLim <- function(x){
    f <- x[is.finite(x)]
    if(length(f))range(f) else range(x)
  }
  h.pixels <- (length(unique(peak.box.counts$chrom))+5)*15
  if(requireNamespace("animint2")){
    viz <- list(
      genome=ggplot2Animint::ggplot()+
        ggplot2Animint::theme_bw()+
        animint2::theme_animint(height=h.pixels)+
        ggplot2Animint::scale_fill_continuous(
          low="grey90", high="red")+
        ggplot2Animint::geom_tile(ggplot2Animint::aes(
          mid.peakEnd/1e6, chrom.fac,
          tooltip=paste(
            count, "peaks in", chrom.box),
          fill=log10(count)),
          clickSelects="chrom.box",
          data=peak.box.counts)+
        ggplot2Animint::scale_x_continuous(
          "position on chromosome (mega bases)")+
        ggplot2Animint::scale_y_discrete("chromosome"),
      peakSize=ggplot2Animint::ggplot()+
        ggplot2Animint::scale_fill_continuous(
          low="grey90", high="blue", na.value="white",
          limits=getLim(log10(peak.box.counts$some.input)))+
        ggplot2Animint::xlab("")+
        ggplot2Animint::geom_point(ggplot2Animint::aes(
          xval, log10(n.samples),
          key=peak.name,
          tooltip=paste(
            n.samples, "samples with a",
            peakBases, "bp peak on",
            peak.name,
            ifelse(
              n.Input==0, "",
              paste("including", n.Input, "Input")),
            sample.counts),
          fill=log10(n.Input),
          href=sprintf(
            "http://genome.ucsc.edu/cgi-bin/hgTracks?position=%s:%d-%d",
            chrom,
            as.integer(peakStart-peakBases*zoom.factor),
            as.integer(min.peakEnd+peakBases*zoom.factor))
        ),
        size=4,
        showSelected="chrom.box",
        chunk_vars="chrom.box",
        data=peak.boxes[, rbind(
          data.table(peak.boxes, x="log10(bases)", xval=log10(peakBases)),
          data.table(
            peak.boxes,
            x="relative position (mega bases)",
            xval=(peakStart-mid.peakEnd)/1e6))])+
        ggplot2Animint::theme_bw()+
        ggplot2Animint::facet_grid(. ~ x, scales="free")+
        ggplot2Animint::theme(panel.margin=grid::unit(0, "lines"))+
        animint2::theme_animint(width=600, height=h.pixels))
    animint2::animint2dir(viz, file.path(set.dir, "figure-genome"))
  }
  figure.png.vec <- Sys.glob(file.path(
    set.dir, "problems", "*", "chunks", "*", "figure-predictions-zoomout.png"))
  if(0 == length(figure.png.vec)){
    chunk.info <- data.table()
    chunks.html <- ""
  }else{
    relative.vec <- sub("/", "", sub(set.dir, "", figure.png.vec))
    zoomin.png.vec <- sub("-zoomout", "", relative.vec)
    chunk.dir.vec <- dirname(zoomin.png.vec)
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
      df <- namedCapture::str_match_named(
        problem,
        g.pos.pattern,
        list(
          chromStart=as.integer,
          chromEnd=as.integer))
      data.frame(problem, df)
    }
    chunk.info <- data.table(
      separate=pos2df(separate.prob.dir.vec),
      chunk=pos2df(chunk.dir.vec),
      zoomin.png=zoomin.png.vec)
    chunk.info[, image := sprintf('
<a href="%s">
  <img src="%s" />
</a>
', zoomin.png.vec, sub(".png$", "-thumb.png", zoomin.png.vec))]
    chunk.info[, chunk := sprintf({
      '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?position=%s">%s</a>'
    }, chunk.problem, chunk.problem)]
    chunk.counts <- chunk.info[, list(chunks=.N), by=separate.problem]
    problems[, labeled.chunks := 0L]
    setkey(chunk.counts, separate.problem)
    setkey(problems, separate.problem)
    problems[chunk.counts, labeled.chunks := chunk.counts$chunks]
    chunks.xt <- xtable(chunk.info[, .(chunk, image)])
    chunks.html <- print(chunks.xt, type="html", sanitize.text.function=identity)
  }
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
  ), by=separate.problem]
  problems[, labels := 0L]
  setkey(label.problem.counts, separate.problem)
  setkey(problems, separate.problem)
  problems[label.problem.counts, labels := label.problem.counts$labels]
  problems[, labeled.regions := 0L]
  problems[label.problem.counts, labeled.regions := label.problem.counts$labeled.regions]
  problems[, problem := sprintf({
    '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?position=%s">%s</a>'
  }, separate.problem, separate.problem)]
  ## add peaks, peak.samples, samples
  separate.counts <- input.pred[, list(
    peak.regions=.N,
    peak.samples=sum(n.samples)
  ), by=list(separate.problem)]
  setkey(separate.counts, separate.problem)
  problems.out <- separate.counts[problems]
  problems.xt <- xtable(problems.out[, list(
    problem,
    peak.regions, peak.samples,
    labeled.chunks, labeled.regions, labels)])
  problems.html <- print(
    problems.xt, type="html", sanitize.text.function=identity)
  html.vec <- c(
    '<title>PeakSegFPOP + PeakSegJoint predictions</title>',
    '<h2>Output: predicted peaks</h2>',
    '<ul>',
    sprintf('<li>
%d genomic regions were found to have
at least one sample with a peak
(<a href="hub.txt">track hub</a>,
 <a href="peaks_summary.tsv">peaks_summary.tsv</a>).
</li>', nrow(input.pred)),
sprintf('<li>
%d peaks were detected overall, across %d samples
(<a href="peaks_matrix.tsv">peaks_matrix.tsv</a>,
<a href="peaks_matrix_height.tsv">peaks_matrix_height.tsv</a>).
</li>',
nrow(joint.peaks.dt),
length(sample.path.vec)),
'<li>
<a href="figure-genome/index.html">Interactive
heatmap and scatterplot of predictions</a>.
</li>',
'</ul>',
'<h2>Input: labeled genomic windows</h2>',
sprintf(
  '%d labeled genomic windows (chunks) were defined in <a href="labels">labels/*.txt</a>',
  nrow(chunk.info)),
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
  peak.mat.dt <- data.table(
    peak=colnames(peak.mat),
    ifelse(t(peak.mat), 1, 0))
  fwrite(
    peak.mat.dt,
    file.path(set.dir, "peaks_matrix.tsv"),
    sep="\t")
  fwrite(
    input.pred,
    file.path(set.dir, "peaks_summary.tsv"),
    sep="\t")
  ## summary
  peaks.bed <- file.path(set.dir, "peaks_summary.bed")
  bed.dt <- input.pred[specificity != "non-specific",]
  max.samples <- max(bed.dt$n.samples)
  bed.dt[, score := as.integer((n.samples/max.samples)*1000) ]
  fwrite(
    bed.dt[, .(
      chrom, peakStart, peakEnd,
      substr(sample.counts, 1, 255),#max characters in bigBed name.
      score)],
    peaks.bed,
    sep="\t",
    col.names=FALSE)
  ## finally, create track hub.
  system(paste("bash", file.path(set.dir, "hub.sh")))
}
