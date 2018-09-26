problem.train <- function
### Train a penalty function that predicts the number of peaks for
### each separate problem. Run this step after computing target
### intervals for each labeled separate problem (problem.target), and
### before peak prediction (problem.predict).
(data.dir.str){
  status <- too.lo <- too.hi <- penalty <- bases <- chromEnd <- chromStart <-
    upper.lim <- lower.lim <- upper.bases <- lower.bases <- ..density.. <-
      prob <- hjust <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  data.dir <- normalizePath(data.dir.str, mustWork=TRUE)
  samples.dir <- file.path(data.dir, "samples")
  model.RData <- file.path(data.dir, "model.RData")
  glob.str <- file.path(
    samples.dir, "*", "*", "problems", "*", "target.tsv")
  cat("Searching for", glob.str, "files for training.\n")
  target.tsv.vec <- Sys.glob(glob.str)
  cat("Found", length(target.tsv.vec), "target.tsv files for training.\n")
  features.list <- list()
  targets.list <- list()
  for(target.tsv.i in seq_along(target.tsv.vec)){
    target.tsv <- target.tsv.vec[[target.tsv.i]]
    problem.dir <- dirname(target.tsv)
    features.tsv <- file.path(problem.dir, "features.tsv")
    if(!file.exists(features.tsv)){
      cat(sprintf("%4d / %4d Computing %s\n", target.tsv.i, length(target.tsv.vec), features.tsv))
      problem.features(problem.dir)
    }
    target.vec <- scan(target.tsv, quiet=TRUE)
    if(any(is.finite(target.vec))){
      features.list[[problem.dir]] <- fread(features.tsv)
      targets.list[[problem.dir]] <- target.vec
    }
  }
  features <- as.matrix(do.call(rbind, features.list))
  targets <- do.call(rbind, targets.list)
  set.seed(1)
  model <- if(nrow(features) < 10){
    some.features <- features[, c("log.quartile.100%", "log.data")]
    cat("Feature matrix:\n")
    print(some.features)
    cat("Target matrix:\n")
    print(unname(targets))
    penaltyLearning::IntervalRegressionUnregularized(
      some.features, targets)
  }else{
    penaltyLearning::IntervalRegressionCV(
      features, targets, verbose=0,
      initial.regularization=1e-4,
      min.observations=nrow(features),
      reg.type=ifelse(nrow(features) < 20, "1sd", "min"))
  }
  model$train.mean.vec <- colMeans(features)
  cat("Learned regularization parameter and weights:\n")
  print(model$pred.param.mat)
  pred.log.penalty <- as.numeric(model$predict(features))
  pred.dt <- data.table(
    problem.dir=rownames(targets),
    too.lo=as.logical(pred.log.penalty < targets[,1]),
    lower.limit=targets[,1],
    pred.log.penalty,
    upper.limit=targets[,2],
    too.hi=as.logical(targets[,2] < pred.log.penalty))
  pred.dt[, status := ifelse(
    too.lo, "low",
                      ifelse(too.hi, "high", "correct"))]
  correct.targets <- pred.dt[status=="correct", ]
  correct.peaks <- correct.targets[!grepl("Input", problem.dir), {
    target_models.tsv <- file.path(problem.dir, "target_models.tsv")
    target.models <- fread(target_models.tsv)
    closest <- target.models[which.min(abs(log(penalty)-pred.log.penalty)),]
    coverage.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
    segments.bed <- paste0(
      coverage.bedGraph, "_penalty=", closest$penalty, "_segments.bed")
    segs <- fread(segments.bed)
    setnames(segs, c("chrom", "chromStart", "chromEnd", "status", "mean"))
    segs[status=="peak", ]
  }, by=problem.dir]
  correct.peaks[, bases := chromEnd-chromStart]
  correct.peaks[, log10.bases := log10(bases)]
  size.model <- correct.peaks[, list(
    mean=mean(log10.bases),
    sd=sd(log10.bases)
  )]
  ## The size.model is constructed by fitting a normal distribution to
  ## the log10(bases) values for all the peaks in the correctly
  ## predicted target intervals in the training data. PeakSegFPOP
  ## typically gives some peaks which are much larger or smaller than
  ## the mean, and these may cause problems for the peak clustering
  ## step. So we use the size.model to exclude peaks which are much
  ## larger or smaller than the mean. The times parameter is the number
  ## of standard deviations past log10(mean) which are allowed. Larger
  ## values of times mean that more peaks will be included in the
  ## separate peak prediction step (times=Inf means that no peaks will
  ## be excluded). Smaller values of times mean that fewer peaks will be
  ## included in the separate peak prediction step (which increases the
  ## risk of false negatives).
  times <- 1
  size.model[, upper.lim := mean + times*sd]
  size.model[, lower.lim := mean - times*sd]
  size.model[, upper.bases := 10^(upper.lim)]
  size.model[, lower.bases := 10^(lower.lim)]
  cat("Train errors:\n")
  print(pred.dt[, list(targets=.N), by=status])
  ## To check if we are extrapolating when predictions are made later,
  ## we save the range of the data 
  model$train.feature.ranges <- apply(
    features[, model$pred.feature.names, drop=FALSE], 2, range)
  ## Plot the size model and limits.
  log10.bases.grid <- correct.peaks[, seq(
    min(log10.bases), max(log10.bases), l=100)]
  normal.dens <- data.table(
    log10.bases=log10.bases.grid,
    prob=size.model[, dnorm(log10.bases.grid, mean, sd)])
  base.labels <- size.model[, {
    log10.bases <- c(lower.lim, mean, upper.lim)
    data.table(
      log10.bases,
      hjust=c(1, 0.5, 0),
      label=scales::comma(round(10^log10.bases)))
  }]
  if(FALSE){
    size.plot <- ggplot()+
      geom_histogram(
        aes(
          log10.bases, ..density..),
        data=correct.peaks)+
      geom_vline(
        aes(xintercept=mean),
        data=size.model,
        size=1, color="red")+
      penaltyLearning::geom_tallrect(
        aes(xmin=lower.lim, xmax=upper.lim),
        data=size.model, fill="red")+
      geom_line(
        aes(log10.bases, prob),
        data=normal.dens, color="red", size=1)+
      geom_text(
        aes(log10.bases, 0, label=label, hjust=hjust),
        data=base.labels, vjust=1)
  }
  cat("Writing model to", model.RData, "\n")
  save(
    model, features, targets,
    size.model,
    correct.peaks,
    file=model.RData)
}

PeakSegFPOP_disk <- structure(function # PeakSegFPOP on disk
### Run the PeakSeg Functional Pruning Optimal Partitioning algorithm,
### using a file on disk (rather than in memory as in
### PeakSegOptimal::PeakSegFPOP) to store the O(N) function piece lists,
### each of size O(log N).
### This is a low-level function that just runs the algo
### and produces the result files (without reading them into R),
### so normal users are recommended to instead use problem.PeakSegFPOP,
### which calls this function then reads the result files into R.
### Finds the optimal change-points using the Poisson loss and the
### PeakSeg constraint. For N data points, the functional pruning
### algorithm is O(N log N) time and disk space, and O(log N) memory.
### It computes the exact
### solution to the following optimization problem. Let Z be an
### N-vector of count data, typically the coverage, number of aligned
### DNA sequence reads in part of the genome
### (the fourth column of bedGraph.file, non-negative integers). Let W
### be an N-vector of positive weights
### (number of bases with the given amount of count/coverage,
### chromEnd - chromStart,
### third column of bedGraph.file - second column). Let penalty
### be a non-negative real number
### (larger for fewer peaks, smaller for more peaks).
### Find the N-vector M of real numbers
### (segment means) and (N-1)-vector C of change-point indicators in
### {-1,0,1} which minimize the penalized Poisson Loss,
### penalty*sum_{i=1}^{N_1} I(c_i=1) + sum_{i=1}^N
### w_i*[m_i-z_i*log(m_i)], subject to constraints: (1) the first
### change is up and the next change is down, etc (sum_{i=1}^t c_i in
### {0,1} for all t<N-1), and (2) the last change is down
### 0=sum_{i=1}^{N-1}c_i, and (3) Every zero-valued change-point
### variable has an equal segment mean after: c_i=0 implies
### m_i=m_{i+1}, (4) every positive-valued change-point variable may
### have an up change after: c_i=1 implies m_i<=m_{i+1}, (5) every
### negative-valued change-point variable may have a down change
### after: c_i=-1 implies m_i>=m_{i+1}. Note that when the equality
### constraints are active for non-zero change-point variables, the
### recovered model is not feasible for the strict inequality
### constraints of the PeakSeg problem, and the optimum of the PeakSeg
### problem is undefined.
(bedGraph.file,
### character scalar: tab-delimited tabular text file with four
### columns: chrom, chromStart, chromEnd, coverage. The algorithm
### creates a large temporary file in the same directory, so make sure
### that there is disk space available on that device.
  pen.str,
### character scalar that can be converted to a numeric scalar via
### as.numeric: non-negative penalty. More penalty means fewer
### peaks. 0 and Inf are OK. Character is required rather than
### numeric, so that the user can reliably find the results in the
### output files, which are in the same directory as bedGraph.file,
### and named using the penalty value,
### e.g. coverage.bedGraph_penalty=136500650856.439_loss.tsv
  allow.free.changes=FALSE
### Allow free changes? If FALSE then the original model with two
### states (up/down) and two changes (up->down is a non-increasing
### change with no penalty, down->up is a non-decreasing change that
### costs penalty). If TRUE then there are two more possible changes:
### up->up non-decreasing for free, down->down non-increasing for
### free.
){
  if(!(
    is.logical(allow.free.changes) &&
    length(allow.free.changes)==1 &&
    is.finite(allow.free.changes)
  )){
    stop("allow.free.changes must be TRUE or FALSE")
  }
  if(!(
    is.character(bedGraph.file) &&
    length(bedGraph.file)==1 &&
    file.exists(bedGraph.file)
  )){
    stop("bedGraph.file must be the name of a data file to segment")
  }
  if(!is.character(pen.str)){
    stop(paste(
      "pen.str must be a character string",
      "that can be converted to a non-negative numeric scalar"
    ))
  }
  penalty <- as.numeric(pen.str)
  if(!(
    is.numeric(penalty) &&
    length(penalty)==1 &&
    0 <= penalty && penalty <= Inf
  )){
    stop("penalty=", penalty, " but it must be a non-negative numeric scalar")
  }
  norm.file <- normalizePath(bedGraph.file, mustWork=TRUE)
  result <- .C(
    "PeakSegFPOP_interface",
    bedGraph.file=as.character(norm.file),
    penalty=pen.str,
    allow.free.changes=allow.free.changes,
    PACKAGE="PeakSegPipeline")
  prefix <- paste0(
    bedGraph.file, "_",
    ifelse(allow.free.changes, "free", ""),
    "penalty=", pen.str)
  result$segments <- paste0(prefix, "_segments.bed")
  result$db <- paste0(prefix, ".db")
  result$loss <- paste0(prefix, "_loss.tsv")
  if(file.size(result$loss)==0){
    stop(
      "unable to write to loss output file ",
      result$loss,
      " (disk is probably full)"
    )
  }
  result
### A list of input parameters (bedGraph.file, penalty) and result
### files (segments, db, loss). 
}, ex=function(){
  
  library(PeakSegPipeline)
  r <- function(chrom, chromStart, chromEnd, coverage){
    data.frame(chrom, chromStart, chromEnd, coverage)
  }
  four <- rbind(
    r("chr1", 0, 10,  2),
    r("chr1", 10, 20, 10),
    r("chr1", 20, 30, 14),
    r("chr1", 30, 40, 13))
  write.table(
    four, tmp <- tempfile(),
    sep="\t", row.names=FALSE, col.names=FALSE)
  names.list <- PeakSegFPOP_disk(tmp, "10.5")
  unlink(names.list$db)
  seg.df <- read.table(names.list$segments)
  names(seg.df) <- c("chrom", "chromStart", "chromEnd", "status", "mean")
  seg.df
  loss.df <- read.table(names.list$loss)
  names(loss.df) <- c(
    "penalty", "segments", "peaks", "bases",
    "mean.pen.cost", "total.cost", "equality.constraints",
    "mean.intervals", "max.intervals")
  loss.df
  
  names.list <- PeakSegFPOP_disk(tmp, "10.5", allow.free.changes=TRUE)
  unlink(names.list$db)
  seg.df <- read.table(names.list$segments)
  names(seg.df) <- c("chrom", "chromStart", "chromEnd", "status", "mean")
  seg.df
  loss.df <- read.table(names.list$loss)
  names(loss.df) <- c(
    "penalty", "segments", "peaks", "bases",
    "mean.pen.cost", "total.cost", "equality.constraints",
    "mean.intervals", "max.intervals")
  loss.df
  
})

problem.predict.allSamples <- function
### Predict for all samples in parallel.
(prob.dir
### project/problems/problemID directory.
 ){
  probs.dir <- dirname(prob.dir)
  set.dir <- dirname(probs.dir)
  problem.name <- basename(prob.dir)
  problem.vec <- Sys.glob(file.path(
    set.dir, "samples", "*", "*", "problems", problem.name))
  peaks.list <- mclapply.or.stop(problem.vec, problem.predict)
  do.call(rbind, peaks.list)
### data.table of predicted peaks.
}

problem.coverage <- function
### Ensure that coverage.bedGraph has been correctly computed for a
### particular genomic segmentation problem.
(problem.dir
### Path to a directory like sampleID/problems/problemID where
### sampleID/coverage.bigWig contains counts of aligned reads in the
### entire genome, and problemID is a chrom range string like
### chr6_dbb_hap3:3491790-3736386 that indicates the genomic
### coordinates of a particular segmentation problem. If
### problemID/coverage.bedGraph does not exist, or its first/last
### lines do not match the expected problemID, then we recreate it
### from sampleID/coverage.bigWig.
){
  chrom <- problemStart <- problemEnd <- count.num.str <- coverage <-
    count.int <- count.int.str <- chromStart <- chromEnd <- J <-
      count <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  problems.dir <- dirname(problem.dir)
  sample.dir <- dirname(problems.dir)
  ## Get problem from directory name.
  pattern <- paste0(
    "(?<chrom>chr[^:]+)",
    ":",
    "(?<problemStart>[0-9]+)",
    "-",
    "(?<problemEnd>[0-9]+)")
  problem.base <- basename(problem.dir)
  problem <- data.table(str_match_named(problem.base, pattern, list(
    problemStart=as.integer,
    problemEnd=as.integer)))
  ## First check if problem/coverage.bedGraph has been created.
  prob.cov.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
  coverage.ok <- tryCatch({
    head.cmd <- paste("head -1", prob.cov.bedGraph) 
    first.cov <- fread(head.cmd)
    setnames(first.cov, c("chrom", "chromStart", "chromEnd", "coverage"))
    tail.cmd <- paste("tail -1", prob.cov.bedGraph)
    last.cov <- fread(tail.cmd)
    setnames(last.cov, c("chrom", "chromStart", "chromEnd", "coverage"))
    last.cov$chromEnd == problem$problemEnd &&
      first.cov$chromStart == problem$problemStart
  }, error=function(e){
    FALSE
  })
  ## If problemID/coverage.bedGraph has already been computed, than we
  ## have nothing to do.
  if(!coverage.ok){
    ## Create problemID/coverage.bedGraph from
    ## sampleID/coverage.bigWig.
    coverage.bigWig <- file.path(sample.dir, "coverage.bigWig")
    cov.cmd <- if(file.exists(coverage.bigWig)){
      problem[, sprintf(
        "bigWigToBedGraph -chrom=%s -start=%d -end=%d %s %s",
        chrom, problemStart, problemEnd,
        coverage.bigWig, prob.cov.bedGraph)]
    }else{
      stop("To compute ", prob.cov.bedGraph,
           " need ", coverage.bigWig,
           " which does not exist.")
    }
    cat(cov.cmd, "\n")
    status <- system(cov.cmd)
    if(status != 0){
      stop("non-zero status code ", status)
    }
    prob.cov <- fread(prob.cov.bedGraph)
    setnames(prob.cov, c("chrom", "chromStart", "chromEnd", "coverage"))
    if(any(prob.cov$coverage < 0)){
      stop("negative coverage in ", prob.cov.bedGraph)
    }
    prob.cov[, count.num.str := paste(coverage)]
    prob.cov[, count.int := as.integer(round(coverage))]
    prob.cov[, count.int.str := paste(count.int)]
    not.int <- prob.cov[count.int.str != count.num.str, ]
    if(nrow(not.int)){
      print(not.int)
      stop("non-integer data in ", prob.cov.bedGraph)
    }
    u.pos <- prob.cov[, sort(unique(c(chromStart, chromEnd)))]
    zero.cov <- data.table(
      chrom=prob.cov$chrom[1],
      chromStart=u.pos[-length(u.pos)],
      chromEnd=u.pos[-1],
      count=0L)
    setkey(zero.cov, chromStart)
    zero.cov[J(prob.cov$chromStart), count := prob.cov$count.int]
    ## The chromStart on the last line of the coverage file should
    ## match the problemEnd, for caching purposes.
    last.end <- zero.cov[.N, chromEnd]
    first.start <- zero.cov[1, chromStart]
    dup.cov <- rbind(if(problem$problemStart==first.start){
      NULL
    }else{
      data.table(
        chrom=prob.cov$chrom[1],
        chromStart=problem$problemStart,
        chromEnd=first.start,
        count=0L)
    }, zero.cov, if(last.end == problem$problemEnd){
      NULL
    }else{
      data.table(
        chrom=prob.cov$chrom[1],
        chromStart=last.end,
        chromEnd=problem$problemEnd,
        count=0L)
    })
    ## dup.cov has chromStart and End the same as problemStart and
    ## end, but maybe has some rows which could be compressed.
    out.cov <- dup.cov[c(diff(count), Inf)!=0]
    out.cov[, chromStart := c(problem$problemStart, chromEnd[-.N])]
    fwrite(
      out.cov,
      prob.cov.bedGraph,
      quote=FALSE,
      sep="\t",
      col.names=FALSE)
  }
### Nothing. If necessary, the bigWigToBedGraph command line program
### is used to create problemID/coverage.bedGraph and then we (1) stop
### if there are any negative or non-integer data and (2) add lines
### with zero counts for missing data.
}

problem.PeakSegFPOP <- structure(function
### Run PeakSegFPOP_disk on one genomic segmentation problem
### directory, and read the result files into R. Actually, this
### function will first check if the result files are already present
### (and consistent), and if so, it will simply read them into R
### (without running PeakSegFPOP_disk) -- this is a caching mechanism
### that can save a lot of time.
(problem.dir,
### Path to a directory like sampleID/problems/problemID which
### contains a coverage.bedGraph file with the aligned read counts for
### one genomic segmentation problem.
 penalty.str,
### Penalty parameter to pass to PeakSegFPOP_disk.
  allow.free.changes=FALSE
### Same as PeakSegFPOP_disk
 ){
  if(!(
    is.logical(allow.free.changes) &&
    length(allow.free.changes)==1 &&
    is.finite(allow.free.changes)
  )){
    stop("allow.free.changes must be TRUE or FALSE")
  }
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  stopifnot(is.character(penalty.str))
  stopifnot(length(penalty.str)==1)
  prob.cov.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
  pre <- paste0(
    prob.cov.bedGraph, "_",
    ifelse(allow.free.changes, "free", ""),
    "penalty=", penalty.str)
  penalty_segments.bed <- paste0(pre, "_segments.bed")
  penalty_loss.tsv <- paste0(pre, "_loss.tsv")
  penalty_timing.tsv <- paste0(pre, "_timing.tsv")
  already.computed <- tryCatch({
    timing <- fread(penalty_timing.tsv)
    setnames(timing, c("penalty", "megabytes", "seconds"))
    first.line <- fread(paste("head -1", penalty_segments.bed))
    setnames(first.line, c("chrom", "chromStart", "chromEnd", "status", "mean"))
    last.line <- fread(paste("tail -1", penalty_segments.bed))
    setnames(last.line, c("chrom", "chromStart", "chromEnd", "status", "mean"))
    penalty.loss <- fread(penalty_loss.tsv)
    setnames(penalty.loss, c(
      "penalty", "segments", "peaks", "bases",
      "mean.pen.cost", "total.cost", "equality.constraints",
      "mean.intervals", "max.intervals"))
    loss.segments.consistent <-
      first.line$chromEnd-last.line$chromStart == penalty.loss$bases
    pattern <- paste0(
      "(?<chrom>chr[^:]+)",
      ":",
      "(?<chromStart>[0-9]+)",
      "-",
      "(?<chromEnd>[0-9]+)")
    problem.base <- basename(problem.dir)
    problem <- str_match_named(problem.base, pattern, list(
      chromStart=as.integer,
      chromEnd=as.integer))
    if(is.na(problem$chrom)){
      stop("problem.dir=", problem.base, " does not match regex ", pattern)
    }
    start.ok <- problem$chromStart == last.line$chromStart
    end.ok <- problem$chromEnd == first.line$chromEnd
    loss.segments.consistent && start.ok && end.ok
  }, error=function(e){
    FALSE
  })
  if(!already.computed){
    penalty.db <- paste0(pre, ".db")
    unlink(penalty.db)#in case interrupted previously.
    seconds <- system.time({
      PeakSegFPOP_disk(prob.cov.bedGraph, penalty.str, allow.free.changes)
    })[["elapsed"]]
    megabytes <- if(file.exists(penalty.db)){
      file.size(penalty.db)/1024/1024
    }else{
      0
    }
    timing <- data.table(
      penalty=as.numeric(penalty.str),
      megabytes,
      seconds)
    write.table(
      timing,
      penalty_timing.tsv,
      row.names=FALSE, col.names=FALSE,
      quote=FALSE, sep="\t")
    unlink(penalty.db)
    penalty.loss <- fread(penalty_loss.tsv)
    setnames(penalty.loss, c(
      "penalty", "segments", "peaks", "bases",
      "mean.pen.cost", "total.cost", "equality.constraints",
      "mean.intervals", "max.intervals"))
  }
  penalty.segs <- fread(penalty_segments.bed)
  setnames(penalty.segs, c("chrom","chromStart", "chromEnd", "status", "mean"))
  sorted.segs <- penalty.segs[order(chromStart)]
  sorted.segs[, diff.status := c(0, diff(status=="peak"))]
  sorted.segs[, peak.i := cumsum(diff.status>0)]
  sorted.segs[, diff.after := c(diff(mean), NA)]
  sorted.segs[, annotation := ifelse(status=="peak", "peakStart", "peakEnd")]
  state.dt <- sorted.segs[, {
    extreme <-.SD[which.max(chromEnd)]
    list(
      stateStart=min(chromStart),
      stateEnd=max(chromEnd),
      extremeMid=extreme[, as.integer((chromEnd+chromStart)/2)],
      extremeMean=extreme$mean
    )}, by=list(peak.i, annotation)]
  penalty.loss[, peaks := nrow(state.dt[annotation=="peakStart"])]
  list(
    states=state.dt,
    segments=sorted.segs,
    loss=penalty.loss,
    timing=timing)
### List of data.tables: segments has one row for every segment in the
### optimal model, timing is one row with the time and disk usage, and
### loss has one row and contains the following columns. penalty=same
### as input, segments=number of segments in optimal model,
### peaks=number of peaks in optimal model, bases=number of positions
### described in bedGraph file, total.cost=total Poisson loss=sum_i
### m_i-z_i*log(m_i)=mean.pen.cost*bases-penalty*peaks,
### mean.pen.cost=mean penalized
### cost=(total.cost+penalty*peaks)/bases, equality.constraints=number
### of adjacent segment means that have equal values in the optimal
### solution, mean.intervals=mean number of intervals/candidate
### changepoints stored in optimal cost functions -- useful for
### characterizing the computational complexity of the algorithm,
### max.intervals=maximum number of intervals.
}, ex=function(){

  library(PeakSegPipeline)
  data(Mono27ac, envir=environment())
  data.dir <- file.path(
    tempfile(),
    "H3K27ac-H3K4me3_TDHAM_BP",
    "samples",
    "Mono1_H3K27ac",
    "S001YW_NCMLS",
    "problems",
    "chr11:60000-580000")
  dir.create(data.dir, recursive=TRUE, showWarnings=FALSE)
  write.table(
    Mono27ac$coverage, file.path(data.dir, "coverage.bedGraph"),
    col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

  ## Compute one model 
  fit <- problem.PeakSegFPOP(data.dir, "5000", allow.free.changes=TRUE)
  ## Visualize that model.
  ann.colors <- c(
    noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")
  library(ggplot2)
  lab.min <- Mono27ac$labels[1, chromStart]
  lab.max <- Mono27ac$labels[.N, chromEnd]
  in.labels <- function(dt){
    dt[lab.min < chromEnd & chromStart < lab.max]
  }

  ggplot()+
    theme_bw()+
    penaltyLearning::geom_tallrect(aes(
      xmin=chromStart, xmax=chromEnd,
      fill=annotation),
      color="grey",
      data=Mono27ac$labels)+
    scale_fill_manual("label", values=ann.colors)+
    geom_step(aes(
      chromStart, count),
      color="grey50",
      data=in.labels(Mono27ac$coverage))+
    geom_segment(aes(
      chromStart, mean,
      xend=chromEnd, yend=mean),
      color="green",
      size=1,
      data=in.labels(fit$segments))+
    geom_segment(aes(
      chromStart, mean,
      xend=chromEnd, yend=mean),
      color="green",
      size=1,
      data=in.labels(fit$segments))+
    ## geom_vline(aes(
    ##   xintercept=chromEnd, linetype=constraint),
    ##   color="green",
    ##   data=in.labels(changes))+
    scale_linetype_manual(values=c(inequality="dotted", equality="solid"))+
    coord_cartesian(xlim=c(lab.min, lab.max))

  peak.y <- -1
  peak.yy <- -2
  sorted.segs <- fit$segments[order(chromStart)]
  sorted.segs[, diff.status := c(0, diff(status=="peak"))]
  sorted.segs[, peak.i := cumsum(diff.status>0)]
  sorted.segs[, diff.after := c(diff(mean), NA)]
  peak.dt <- sorted.segs[0 < peak.i, {
    is.peak <- status=="peak"
    list(
      chromStart=min(chromStart[is.peak]),
      chromEnd=max(chromEnd[is.peak]),
      max.diff.after=chromEnd[which.max(diff.after)],
      min.diff.after=chromEnd[which.min(diff.after)],
      last.peak.mean=mean[max(which(is.peak))]
    )}, by=list(peak.i)]
  gg <- ggplot()+
    theme_bw()+
    penaltyLearning::geom_tallrect(aes(
      xmin=chromStart, xmax=chromEnd,
      fill=annotation),
      color="grey",
      data=Mono27ac$labels)+
    scale_fill_manual("label", values=ann.colors)+
    geom_step(aes(
      chromStart, count),
      color="grey50",
      data=Mono27ac$coverage)+
    geom_segment(aes(
      chromStart, mean,
      xend=chromEnd, yend=mean),
      color="green",
      size=1,
      data=sorted.segs)+
    geom_segment(aes(
      chromStart, mean,
      xend=chromEnd, yend=mean),
      color="green",
      size=1,
      data=sorted.segs)+
    geom_segment(aes(
      chromStart, peak.y,
      xend=chromEnd, yend=peak.y),
      color="deepskyblue",
      size=1.5,
      data=peak.dt)+
    geom_label(aes(
      chromStart, peak.y, label=peak.i),
      color="deepskyblue",
      data=peak.dt)+
    geom_segment(aes(
      max.diff.after, peak.yy,
      xend=min.diff.after, yend=peak.yy),
      color="deepskyblue",
      size=1.5,
      data=peak.dt)+
    geom_label(aes(
      max.diff.after, peak.yy, label=peak.i),
      color="deepskyblue",
      data=peak.dt)+
    geom_label(aes(
      chromEnd, last.peak.mean, label=peak.i),
      color="green",
      data=peak.dt)+
    scale_linetype_manual(values=c(inequality="dotted", equality="solid"))
  print(gg)

  gg+
    coord_cartesian(xlim=c(2e5, 3e5))

  gg+
    coord_cartesian(xlim=c(4.5e5, 5.5e5))


  peak.y <- -1
  peak.yy <- -2
  gg <- ggplot()+
    theme_bw()+
    penaltyLearning::geom_tallrect(aes(
      xmin=chromStart, xmax=chromEnd,
      fill=annotation),
      color="grey",
      data=Mono27ac$labels)+
    scale_fill_manual("label", values=ann.colors)+
    geom_step(aes(
      chromStart, count),
      color="grey50",
      data=Mono27ac$coverage)+
    geom_segment(aes(
      chromStart, mean,
      xend=chromEnd, yend=mean),
      color="green",
      size=1,
      data=sorted.segs)+
    geom_segment(aes(
      chromStart, mean,
      xend=chromEnd, yend=mean),
      color="green",
      size=1,
      data=sorted.segs)+
    geom_segment(aes(
      chromStart, peak.y,
      xend=chromEnd, yend=peak.y),
      color="deepskyblue",
      size=1.5,
      data=peak.dt)+
    geom_label(aes(
      chromStart, peak.y, label=peak.i),
      color="deepskyblue",
      data=peak.dt)+
    geom_segment(aes(
      max.diff.after, peak.yy,
      xend=min.diff.after, yend=peak.yy),
      color="deepskyblue",
      size=1.5,
      data=peak.dt)+
    geom_label(aes(
      max.diff.after, peak.yy, label=peak.i),
      color="deepskyblue",
      data=peak.dt)+
    geom_label(aes(
      chromEnd, last.peak.mean, label=peak.i),
      color="green",
      data=peak.dt)+
    scale_linetype_manual(values=c(inequality="dotted", equality="solid"))
  print(gg)

  gg+
    coord_cartesian(xlim=c(2e5, 3e5))

  gg+
    coord_cartesian(xlim=c(4.5e5, 5.5e5))

  
  ## Compute one two change model 
  fit <- problem.PeakSegFPOP(data.dir, "3000", allow.free.changes=FALSE)
  ## Visualize that model.
  ann.colors <- c(
    noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")
  library(ggplot2)
  lab.min <- Mono27ac$labels[1, chromStart]
  lab.max <- Mono27ac$labels[.N, chromEnd]
  in.labels <- function(dt){
    dt[lab.min < chromEnd & chromStart < lab.max]
  }
  changes <- fit$segments[, list(
    constraint=ifelse(diff(mean)==0, "equality", "inequality"),
    chromStart=chromEnd[-1],
    chromEnd=chromEnd[-1])]

  ggplot()+
    theme_bw()+
    penaltyLearning::geom_tallrect(aes(
      xmin=chromStart, xmax=chromEnd,
      fill=annotation),
      color="grey",
      data=Mono27ac$labels)+
    scale_fill_manual("label", values=ann.colors)+
    geom_step(aes(
      chromStart, count),
      color="grey50",
      data=in.labels(Mono27ac$coverage))+
    geom_segment(aes(
      chromStart, mean,
      xend=chromEnd, yend=mean),
      color="green",
      size=1,
      data=in.labels(fit$segments))+
    geom_segment(aes(
      chromStart, mean,
      xend=chromEnd, yend=mean),
      color="green",
      size=1,
      data=in.labels(fit$segments))+
    ## geom_vline(aes(
    ##   xintercept=chromEnd, linetype=constraint),
    ##   color="green",
    ##   data=in.labels(changes))+
    scale_linetype_manual(values=c(inequality="dotted", equality="solid"))+
    coord_cartesian(xlim=c(lab.min, lab.max))

  peak.y <- -1
  peak.dt <- fit$segments[status=="peak"]
  peak.dt[, peak.i := seq_along(chrom)]
  gg <- ggplot()+
    theme_bw()+
    penaltyLearning::geom_tallrect(aes(
      xmin=chromStart, xmax=chromEnd,
      fill=annotation),
      color="grey",
      data=Mono27ac$labels)+
    scale_fill_manual("label", values=ann.colors)+
    geom_step(aes(
      chromStart, count),
      color="grey50",
      data=Mono27ac$coverage)+
    geom_segment(aes(
      chromStart, mean,
      xend=chromEnd, yend=mean),
      color="green",
      size=1,
      data=sorted.segs)+
    geom_segment(aes(
      chromStart, mean,
      xend=chromEnd, yend=mean),
      color="green",
      size=1,
      data=sorted.segs)+
    geom_vline(aes(
      xintercept=chromEnd, linetype=constraint),
      color="green",
      data=changes)+
    geom_segment(aes(
      chromStart, peak.y,
      xend=chromEnd, yend=peak.y),
      color="deepskyblue",
      size=1.5,
      data=peak.dt)+
    geom_label(aes(
      chromStart, peak.y, label=peak.i),
      color="deepskyblue",
      data=peak.dt)+
    scale_linetype_manual(values=c(inequality="dotted", equality="solid"))
  print(gg)

  gg+
    coord_cartesian(xlim=c(2e5, 3e5))

  gg+
    coord_cartesian(xlim=c(4.5e5, 5.5e5))

  
})

problem.sequentialSearch <- structure(function
### Compute the most likely peak model with at most the number of
### peaks given by peaks.int. This function repeated calls
### problem.PeakSegFPOP with different penalty values, until either
### (1) it finds the peaks.int model, or (2) it concludes that there
### is no peaks.int model, in which case it returns the next simplest
### model (with fewer peaks than peaks.int).
(problem.dir,
### problemID directory in which coverage.bedGraph has already been
### computed. If there is a labels.bed file then the number of
### incorrect labels will be computed in order to find the target
### interval of minimal error penalty values.
  peaks.int,
### int: target number of peaks.
  verbose=0,
### Print messages?
  allow.free.changes=FALSE
### TRUE for four edge state graph with free changes, FALSE for two
### edges.
){
  stopifnot(
    is.integer(peaks.int) &&
    length(peaks.int)==1 &&
    0 <= peaks.int)
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  model.list <- list()
  next.pen <- c(0, Inf)
  iteration <- 0
  under <- over <- data.table(peaks=NA)
  while(length(next.pen)){
    if(verbose)cat(
      "Next =", paste(next.pen, collapse=", "),
      "mc.cores=", getOption("mc.cores"),
      "\n")
    next.str <- paste(next.pen)
    iteration <- iteration+1
    model.list[next.str] <- mclapply.or.stop(
      next.str, function(penalty.str){
        L <- problem.PeakSegFPOP(
          problem.dir,
          penalty.str,
          allow.free.changes=allow.free.changes)
        L$loss$seconds <- L$timing$seconds
        L$loss$megabytes <- L$timing$megabytes
        L$loss$iteration <- iteration
        L$loss$under <- under$peaks
        L$loss$over <- over$peaks
        L
      }
    )
    if(iteration==1){
      under <- model.list[["Inf"]]$loss
      over <- model.list[["0"]]$loss
      max.peaks <- floor((over$bases-1)/2)
      if(max.peaks < peaks.int){
        stop(
          "peaks.int=",
          peaks.int,
          " but max=",
          max.peaks,
          " peaks for N=",
          over$bases,
          " data")
      }
    }else{
      Mnew <- model.list[[next.str]]$loss
      if(Mnew$peaks %in% c(under$peaks, over$peaks)){## not a new model.
        candidate <- under ##pick the simpler one.
        next.pen <- NULL
      }else{#new model.
        if(Mnew$peaks < peaks.int){
          under <- Mnew
        }else{
          over <- Mnew
        }
      }
    }
    if(peaks.int==under$peaks){
      candidate <- under
      next.pen <- NULL
    }
    if(peaks.int==over$peaks){
      candidate <- over
      next.pen <- NULL
    }
    if(!is.null(next.pen)){
      next.pen <- (over$total.cost-under$total.cost)/(under$peaks-over$peaks)
      if(next.pen<0){
        ## sometimes happens for a large number of peaks -- cost is
        ## numerically unstable so we don't get a good penalty to try --
        ## anyways these models are way too big, so just return under.
        candidate <- under
        next.pen <- NULL
      }
    }
  }#while(!is.null(pen))
  out <- model.list[[paste(candidate$penalty)]]
  loss.list <- lapply(model.list, "[[", "loss")
  out$others <- do.call(rbind, loss.list)[order(iteration)]
  out
### Same result list from problem.PeakSegFPOP, with an additional
### component "others" describing the other models that were computed
### before finding the optimal model with peaks.int (or fewer)
### peaks. Additional loss columns are as follows: under=number of
### peaks in smaller model during binary search; over=number of peaks
### in larger model during binary search; iteration=number of times
### PeakSegFPOP has been run.
}, ex=function(){

  ## Create simple 6 point data set discussed in supplementary
  ## materials. GFPOP/GPDPA computes up-down model with 2 peaks, but
  ## neither CDPA (PeakSegDP::cDPA) nor PDPA (jointseg)
  r <- function(chrom, chromStart, chromEnd, coverage){
    data.frame(chrom, chromStart, chromEnd, coverage)
  }
  supp <- rbind(
    r("chr1", 0, 1,  3),
    r("chr1", 1, 2, 9),
    r("chr1", 2, 3, 18),
    r("chr1", 3, 4, 15),
    r("chr1", 4, 5, 20),
    r("chr1", 5, 6, 2)
  )
  data.dir <- file.path(tempfile(), "chr1:0-6")
  dir.create(data.dir, recursive=TRUE)
  write.table(
    supp, file.path(data.dir, "coverage.bedGraph"),
    sep="\t", row.names=FALSE, col.names=FALSE)

  ## Compute optimal up-down model with 2 peaks via sequential search.
  fit <- PeakSegPipeline::problem.sequentialSearch(data.dir, 2L)

  library(ggplot2)
  ggplot()+
    theme_bw()+
    geom_point(aes(
      chromEnd, coverage),
      data=supp)+
    geom_segment(aes(
      chromStart+0.5, mean,
      xend=chromEnd+0.5, yend=mean),
      data=fit$segments,
      color="green")
  
})


problem.features <- function
### Compute features for one segmentation problem.
(problem.dir
### problemID directory with problemID/coverage.bedGraph.
 ){
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  coverage <- fread(file.path(problem.dir, "coverage.bedGraph"))
  setnames(coverage, c("chrom", "chromStart", "chromEnd", "count"))
  bases <- with(coverage, chromEnd-chromStart)
  long <- rep(coverage$count, bases)
  diff.vec <- abs(diff(long))
  feature.vec <- c(
    quartile=quantile(long),
    mean=mean(long),
    sd=sd(long),
    bases=sum(bases),
    data=nrow(coverage))
  log.features <- suppressWarnings(c(
    feature.vec,
    `log+1`=log(feature.vec+1),
    log=log(feature.vec),
    log.log=log(log(feature.vec))))
  feature.dt <- data.table(t(log.features))
  write.table(
    feature.dt,
    file.path(problem.dir, "features.tsv"),
    quote=FALSE,
    row.names=FALSE,
    col.names=TRUE,
    sep="\t")
### Nothing, but writes problemID/features.tsv if it does not exist
### already.
}

problem.target <- structure(function
### Compute target interval for a segmentation problem. This function
### repeated calls problem.PeakSegFPOP with different penalty values,
### until it finds an interval of penalty values with minimal label
### error. The calls to PeakSegFPOP are parallelized using mclapply if
### you set options(mc.cores).
(problem.dir,
### problemID directory in which coverage.bedGraph has already been
### computed. If there is a labels.bed file then the number of
### incorrect labels will be computed in order to find the target
### interval of minimal error penalty values.
  minutes.limit=NULL,
### Time limit; the search will stop at a sub-optimal target interval
### if this many minutes has elapsed. Useful for testing environments
### with build time limits (travis). Default NULL means to use the
### value in option PeakSegPipeline.problem.target.minutes (or Inf if
### that option is not set).
  verbose=0
 ){
  status <- peaks <- errors <- fp <- fn <- penalty <- max.log.lambda <-
    min.log.lambda <- penalty <- . <- done <- total.cost <- mean.pen.cost <-
      bases <- no.next <- is.min <- min.err.interval <- max.lambda <-
        already.computed <- is.other <- dist <- min.lambda <- log.size <-
          mid.lambda <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  if(is.null(minutes.limit)){
    ## here rather than in the arguments in order to avoid Rd NOTE
    ## about lines wider than 90 characters.
    minutes.limit <- getOption("PeakSegPipeline.problem.target.minutes", Inf)
  }
  seconds.start <- as.numeric(Sys.time())
  stopifnot(is.numeric(minutes.limit))
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  problem.coverage(problem.dir)
  ## Check if problem/labels.bed exists.
  problem.labels <- tryCatch({
    prob.lab.bed <- file.path(problem.dir, "labels.bed")
    problem.labels <- fread(prob.lab.bed)
    setnames(problem.labels, c("chrom", "chromStart", "chromEnd", "annotation"))
    problem.labels
  }, error=function(e){
    data.frame(
      chrom=character(),
      chromStart=integer(),
      chromEnd=integer(),
      annotation=character())
  })
  ## Compute the label error for one penalty parameter.
  getError <- function(penalty.str){
    stopifnot(is.character(penalty.str))
    stopifnot(length(penalty.str) == 1)
    result <- problem.PeakSegFPOP(problem.dir, penalty.str)
    penalty.peaks <- result$segments[status=="peak",]
    tryCatch({
      penalty.error <- PeakErrorChrom(penalty.peaks, problem.labels)
    }, error=function(e){
      stop("try deleting _segments.bed and recomputing, error computing number of incorrect labels: ", e)
    })
    with(penalty.error, data.table(
      iteration,
      result$loss,
      fn=sum(fn),
      fp=sum(fp)))
  }  
  ## Also compute feature vector here so train is faster later.
  problem.features(problem.dir)
  error.list <- list()
  next.pen <- c(0, Inf)
  iteration <- 0
  last.target.vec <- c(-Inf, Inf)
  target.result.list <- list()
  while(length(next.pen)){
    if(verbose)cat(
      "Next =", paste(next.pen, collapse=", "),
      "mc.cores=", getOption("mc.cores"),
      "\n")
    next.str <- paste(next.pen)
    iteration <- iteration+1
    error.list[next.str] <- mclapply.or.stop(next.str, getError)
    error.dt <- do.call(rbind, error.list)[order(penalty)]
    if(!is.numeric(error.dt$penalty)){
      stop("penalty column is not numeric -- check loss in _loss.tsv files")
    }
    error.dt[, errors := fp+fn]
    if(verbose)print(error.dt[,.(penalty, peaks, status, fp, fn, errors)])
    unique.peaks <- error.dt[, data.table(
      .SD[which.min(iteration)],
      penalties=.N
    ), by=list(peaks)]
    path.dt <- data.table(penaltyLearning::modelSelection(
      unique.peaks, "total.cost", "peaks"))
    path.dt[, next.pen := max.lambda]
    path.dt[, already.computed := next.pen %in% names(error.list)]
    path.dt[, no.next := c(diff(peaks) == -1, NA)]
    path.dt[, done := already.computed | no.next]
    path.dt[, is.min := errors==min(errors)]
    path.dt[, min.err.interval := cumsum(ifelse(
      c(is.min[1], diff(is.min))==1, 1, 0))]
    other.candidates <- path.dt[which(0<diff(fn) & diff(fp)<0)]
    interval.dt <- path.dt[is.min==TRUE, {
      i <- if(1 == .N || 0 == errors[1]){
        ## No middle candidate if there is only one model in the
        ## interval, or if there are no errors.
        NA
      }else{
        d <- data.table(
          i=1:(.N-1),
          ## do not attempt to explore other.candidates -- try
          ## different ones!
          is.other=next.pen[-.N] %in% other.candidates$next.pen,
          dist=diff(max.log.lambda)+diff(min.log.lambda),
          done=done[-.N])
        d[is.other==FALSE & done==FALSE, i[which.max(dist)]]
      }
      if(length(i)==0)i <- NA
      data.table(
        min.lambda=min.lambda[1],
        min.log.lambda=min.log.lambda[1],
        mid.lambda=max.lambda[i],
        max.lambda=max.lambda[.N],
        max.log.lambda=max.log.lambda[.N],
        log.size=max.log.lambda[.N]-min.log.lambda[1]
        )
    }, by=list(min.err.interval)]
    largest.interval <- interval.dt[which.max(log.size)]
    target.vec <- largest.interval[, c(min.log.lambda, max.log.lambda)]
    write(target.vec, file.path(problem.dir, "target.tsv"), sep="\t")
    diff.target.vec <- target.vec-last.target.vec
    last.target.vec <- target.vec
    target.result.list[[paste(iteration)]] <- largest.interval[, data.table(
      iteration,
      min.log.lambda,
      max.log.lambda)]
    target.lambda <- largest.interval[, c(min.lambda, max.lambda)]
    error.candidates <- path.dt[next.pen %in% target.lambda]
    stopping.candidates <- rbind(error.candidates, other.candidates)[done==FALSE]
    seconds.now <- as.numeric(Sys.time())
    minutes.elapsed <- (seconds.now-seconds.start)/60
    if(verbose)cat(sprintf(
      "%f minutes elapsed / %f limit\nTarget interval: %f %f change: %f %f\n",
      minutes.elapsed, minutes.limit,
      target.vec[1], target.vec[2],
      diff.target.vec[1], diff.target.vec[2]))
    next.pen <- if(minutes.elapsed < minutes.limit && nrow(stopping.candidates)){
      lambda.vec <- interval.dt[, c(min.lambda, mid.lambda, max.lambda)]
      interval.candidates <- path.dt[next.pen %in% lambda.vec][done==FALSE]
      unique(rbind(stopping.candidates, interval.candidates)$next.pen)
    }
  }#while(!is.null(pen))
  write.table(
    error.dt,
    file.path(problem.dir, "target_models.tsv"),
    sep="\t",
    quote=FALSE,
    row.names=FALSE,
    col.names=TRUE)
  list(
    target=target.vec,
    target.iterations=do.call(rbind, target.result.list),
    models=error.dt)
### List of info related to target interval computation: target is the
### interval of log(penalty) values that achieve minimum incorrect
### labels (numeric vector of length 2), target.iterations is a
### data.table with target intervals as a function of iteration,
### models is a data.table with one row per model for which the label
### error was computed.
}, ex=function(){

  library(PeakSegPipeline)
  data(Mono27ac, envir=environment())
  ## Write the Mono27ac data set to disk.
  data.dir <- file.path(
    tempfile(),
    "H3K27ac-H3K4me3_TDHAM_BP",
    "samples",
    "Mono1_H3K27ac",
    "S001YW_NCMLS",
    "problems",
    "chr11:60000-580000")
  dir.create(data.dir, recursive=TRUE, showWarnings=FALSE)
  write.table(
    Mono27ac$labels, file.path(data.dir, "labels.bed"),
    col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
  write.table(
    Mono27ac$coverage, file.path(data.dir, "coverage.bedGraph"),
    col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

  ## Compute target interval. Specifying minutes.limit stops the
  ## optimization after that number of minutes, resulting in an
  ## imprecise target interval, but saving time (to avoid NOTE on
  ## CRAN).
  target.list <- problem.target(data.dir, minutes.limit=0.05)

  ## These are all the models computed in order to find the target
  ## interval.
  print(target.list$models[, list(
    penalty, log.penalty=log(penalty), peaks, total.cost, fn, fp, errors)])

  ## This is the target interval in log(penalty) values.
  print(target.list$target)

})

problem.predict <- function
### Predict peaks for a genomic segmentation problem.
(problem.dir
### project/samples/groupID/sampleID/problems/problemID.
 ){
  model <- status <- chromEnd <- chromStart <- size.model <- lower.bases <-
    upper.bases <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  problems.dir <- dirname(problem.dir)
  sample.dir <- dirname(problems.dir)
  group.dir <- dirname(sample.dir)
  samples.dir <- dirname(group.dir)
  data.dir <- dirname(samples.dir)
  sample.id <- basename(sample.dir)
  sample.group <- basename(group.dir)
  cov.result <- try(problem.coverage(problem.dir))
  if(inherits(cov.result, "try-error")){
    cat("Could not compute coverage in", problem.dir,
        "so not predicting peaks.\n")
    return(NULL)
  }
  features.tsv <- file.path(problem.dir, "features.tsv")
  is.computed <- if(file.exists(features.tsv)){
    TRUE
  }else{
    tryCatch({
      problem.features(problem.dir)
      cat(sprintf("Computed %s\n", features.tsv))
      TRUE
    }, error=function(e){
      FALSE
    })
  }
  if(!is.computed){
    cat("Unable to compute", features.tsv, "so not predicting.\n")
    return(NULL)
  }
  features <- fread(features.tsv)
  stopifnot(nrow(features)==1)
  feature.mat <- as.matrix(features)
  model.RData <- file.path(data.dir, "model.RData")
  load(model.RData)
  if(length(feature.mat)==length(model$train.mean.vec)){
    is.bad <- !is.finite(feature.mat)
    feature.mat[is.bad] <- model$train.mean.vec[is.bad]
  }else{
    stop(
      "feature.mat has ", length(feature.mat),
      " columns but ", length(model$train.mean.vec),
      " features were used to train the model")
  }
  pred.penalty <- as.numeric(exp(model$predict(feature.mat)))
  stopifnot(length(pred.penalty)==1)
  stopifnot(is.finite(pred.penalty))
  n.features <- length(model$pred.feature.names)
  cat(paste0(
    "Predicting penalty=", pred.penalty,
    " log(penalty)=", log(pred.penalty),
    " based on ", n.features,
    " feature", ifelse(n.features==1, "", "s"),
    ".\n"))
  pen.str <- paste(pred.penalty)
  result <- problem.PeakSegFPOP(problem.dir, pen.str)
  all.peaks <- result$segments[status=="peak", ]
  bases.vec <- all.peaks[, chromEnd-chromStart]
  in.range <- size.model[, lower.bases < bases.vec & bases.vec < upper.bases]
  peaks <- all.peaks[in.range, ]
  ## save peaks.
  peaks.bed <- file.path(problem.dir, "peaks.bed")
  cat(
    "Writing ", peaks.bed,
    " with ", nrow(peaks),
    " peak", ifelse(nrow(peaks)==1, "", "s"),
      ".\n", sep="")
  write.table(
    peaks,
    peaks.bed,
    quote=FALSE,
    sep="\t",
    col.names=FALSE,
    row.names=FALSE)
  if(nrow(peaks)){
    data.table(sample.id, sample.group, peaks)
  }
### data.table of peak predictions.
}
