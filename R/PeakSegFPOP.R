problem.tempfile <- function
### Create a (problem,penalty)-specific temporary file name to pass to
### PeakSegFPOP_dir as the cost function database.
(problem.dir,
### full path to problem directory.
  pen.str
### penalty string.
){
  problems.dir <- dirname(problem.dir)
  sample.dir <- dirname(problems.dir)
  group.dir <- dirname(sample.dir)
  samples.dir <- dirname(group.dir)
  data.dir <- dirname(samples.dir)
  sample.id <- basename(sample.dir)
  sample.group <- basename(group.dir)
  data.name <- basename(data.dir)
  data.prob.pen <- paste(
    data.name, sample.group, sample.id,
    basename(problem.dir),
    pen.str,
    sep="_")
  file.path(tempdir(), data.prob.pen)
### character: temporary file name.
}

problem.train <- function
### Train a penalty function that predicts the number of peaks for
### each separate problem. Run this step after computing target
### intervals for each labeled separate problem (problem.target), and
### before peak prediction (problem.predict).
(data.dir.str,
  verbose=getOption("PeakSegPipeline.verbose", 1)
){
  status <- too.lo <- too.hi <- penalty <- bases <- chromEnd <- chromStart <-
    upper.lim <- lower.lim <- upper.bases <- lower.bases <- ..density.. <-
      prob <- hjust <- . <- min.log.lambda <- max.log.lambda <- median.bases <-
        cached.bases <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  data.dir <- normalizePath(data.dir.str, mustWork=TRUE)
  samples.dir <- file.path(data.dir, "samples")
  model.RData <- file.path(data.dir, "model.RData")
  train_data.csv <- file.path(data.dir, "train_data.csv")
  glob.str <- file.path(
    samples.dir, "*", "*", "problems", "*", "target.tsv")
  if(verbose)cat(
    "Searching for", glob.str, "files for training.\n")
  target.tsv.vec <- Sys.glob(glob.str)
  if(verbose)cat(
    "Found", length(target.tsv.vec), "new target.tsv files for training.\n")
  train.dt <- if(file.exists(train_data.csv)){
    train.all <- fread(file=train_data.csv)
    train.all[!problem.dir %in% dirname(target.tsv.vec)]
  }else{
    data.table()
  }
  train.dt.list <- list(cached=train.dt)
  for(target.tsv.i in seq_along(target.tsv.vec)){
    target.tsv <- target.tsv.vec[[target.tsv.i]]
    problem.dir <- dirname(target.tsv)
    if(verbose)cat(sprintf(
      "%4d / %4d Reading features/targets %s\n",
      target.tsv.i, length(target.tsv.vec), problem.dir))
    features.tsv <- file.path(problem.dir, "features.tsv")
    if(!file.exists(features.tsv)){
      problem.features(problem.dir)
    }
    feature.row <- fread(file=features.tsv)
    target.row <- fread(
      file=target.tsv,
      col.names=c("min.log.lambda", "max.log.lambda"))
    train.dt.list[[problem.dir]] <- data.table(
      problem.dir, target.row, feature.row)
  }
  train.dt <- do.call(rbind, train.dt.list)
  fwrite(train.dt, train_data.csv)
  unlink(target.tsv.vec)
  some.finite <- train.dt[
    is.finite(min.log.lambda) | is.finite(max.log.lambda)]
  features <- as.matrix(some.finite[, -(1:3), with=FALSE])
  targets <- as.matrix(some.finite[, .(min.log.lambda, max.log.lambda)])
  set.seed(1)
  model <- if(nrow(features) < 10){
    some.features <- features[, c("log.quartile.100%", "log.data")]
    if(verbose)cat("Feature matrix:\n")
    if(verbose)print(some.features)
    if(verbose)cat("Target matrix:\n")
    if(verbose)print(unname(targets))
    penaltyLearning::IntervalRegressionUnregularized(
      some.features, targets)
  }else{
    penaltyLearning::IntervalRegressionCV(
      features, targets, verbose=getOption("PeakSegPipeline.verbose", 1),
      initial.regularization=1e-4,
      min.observations=nrow(features),
      reg.type=ifelse(nrow(features) < 20, "1sd", "min"))
  }
  model$train.mean.vec <- colMeans(features)
  if(verbose)cat("Learned regularization parameter and weights:\n")
  if(verbose)print(model$pred.param.mat)
  pred.log.penalty <- as.numeric(model$predict(features))
  pred.dt <- data.table(
    problem.dir=some.finite$problem.dir,
    too.lo=as.logical(pred.log.penalty < targets[,1]),
    lower.limit=targets[,1],
    pred.log.penalty,
    upper.limit=targets[,2],
    too.hi=as.logical(targets[,2] < pred.log.penalty))
  pred.dt[, status := ifelse(
    too.lo, "low",
    ifelse(too.hi, "high", "correct"))]
  correct.targets <- pred.dt[status=="correct"]
  correct.peak.stats <- correct.targets[!grepl("Input", problem.dir), .(
    problem.dir,
    median.bases=NA_real_,
    pred.log.penalty)]
  correct_peaks.csv <- file.path(data.dir, "correct_peaks.csv")
  if(file.exists(correct_peaks.csv)){
    correct.cache <- fread(file=correct_peaks.csv)
    if(nrow(correct.cache)){
      correct.peak.stats[
        correct.cache,
        median.bases := as.numeric(cached.bases),
        on="problem.dir"]
    }
  }
  correct.peak.stats[is.na(median.bases), median.bases := {
    models.rds <- file.path(problem.dir, "models.rds")
    models.dt <- readRDS(models.rds)
    closest <- models.dt[which.min(abs(log(penalty)-pred.log.penalty))]
    segs.dt <- closest$segments.dt[[1]]
    segs.dt[status=="peak", as.numeric(median(chromEnd-chromStart))]
  }, by=problem.dir]
  new.cache <- correct.peak.stats[, .(problem.dir, cached.bases=median.bases)]
  fwrite(new.cache, correct_peaks.csv)
  correct.peak.stats[, log10.bases := log10(median.bases)]
  size.model <- correct.peak.stats[, list(
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
  if(verbose)cat("Train errors:\n")
  if(verbose)print(pred.dt[, list(targets=.N), by=status])
  ## To check if we are extrapolating when predictions are made later,
  ## we save the range of the data
  model$train.feature.ranges <- apply(
    features[, model$pred.feature.names, drop=FALSE], 2, range)
  ## Plot the size model and limits.
  log10.bases.grid <- correct.peak.stats[, seq(
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
        data=correct.peak.stats)+
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
  if(verbose)cat("Writing model to", model.RData, "\n")
  save(
    model, features, targets,
    size.model,
    file=model.RData)
}

problem.predict.allSamples <- function
### Predict for all samples, parallelized over problems via
### future.apply::future_lapply.
(prob.dir
### project/problems/problemID directory.
 ){
  probs.dir <- dirname(prob.dir)
  set.dir <- dirname(probs.dir)
  problem.name <- basename(prob.dir)
  sample.dir.vec <- Sys.glob(file.path(
    set.dir, "samples", "*", "*"))
  prob.name <- basename(prob.dir)
  problem.vec <- file.path(sample.dir.vec, "problems", prob.name)
  peaks.list <- future.apply::future_lapply(problem.vec, problem.predict)
  do.call(rbind, peaks.list)
### data.table of predicted peaks.
}

problem.table <- function
### Convert a path with a problem string to a data.table.
(problem.dir
### path with a problem string, e.g. chrX:6000-1000000
){
  nc::capture_first_vec(
    problem.dir,
    chrom="chr[^-:]+",
    "[-:]",
    problemStart="[0-9]+", as.integer,
    "-",
    problemEnd="[0-9]+", as.integer)
### data.table with columns chrom, problemStart, problemEnd.
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
  dir.create(problem.dir, showWarnings=FALSE, recursive=TRUE)
  problems.dir <- dirname(problem.dir)
  sample.dir <- dirname(problems.dir)
  problem <- problem.table(problem.dir)
  ## First check if problem/coverage.bedGraph has been created.
  prob.cov.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
  coverage.ok <- tryCatch({
    first.cov <- suppressWarnings(fread.first(prob.cov.bedGraph, c(
      "chrom", "chromStart", "chromEnd", "coverage")))
    last.cov <- suppressWarnings(fread.last(prob.cov.bedGraph, c(
      "chrom", "chromStart", "chromEnd", "coverage")))
    is.integer(first.cov$chromStart) &&
      is.integer(last.cov$chromEnd) &&
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
    if(!file.exists(coverage.bigWig)){
      stop("To compute ", prob.cov.bedGraph,
           " need ", coverage.bigWig,
           " which does not exist.")
    }
    prob.cov <- problem[, readBigWig(
      coverage.bigWig,
      chrom,
      problemStart,
      problemEnd)]
    if(nrow(prob.cov)==0){
      stop(
        "coverage/count data file ",
        prob.cov.bedGraph,
        " is empty; typically this happens when ",
        coverage.bigWig,
        " has no data in this genomic region")
    }
    if(any(prob.cov$count < 0)){
      stop("negative coverage in ", prob.cov.bedGraph)
    }
    prob.cov[, count.num.str := paste(count)]
    prob.cov[, count.int := as.integer(round(count))]
    prob.cov[, count.int.str := paste(count.int)]
    not.int <- prob.cov[count.int.str != count.num.str, ]
    if(nrow(not.int)){
      print(not.int)
      stop("non-integer data in ", prob.cov.bedGraph)
    }
    u.pos <- prob.cov[, sort(unique(c(chromStart, chromEnd)))]
    zero.cov <- data.table(
      chrom=problem$chrom,
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
        chrom=problem$chrom,
        chromStart=problem$problemStart,
        chromEnd=first.start,
        count=0L)
    }, zero.cov, if(last.end == problem$problemEnd){
      NULL
    }else{
      data.table(
        chrom=problem$chrom,
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
  problem
### problem data.table. If necessary, the bigWigToBedGraph command
### line program is used to create problemID/coverage.bedGraph and
### then we (1) stop if there are any negative or non-integer data and
### (2) add lines with zero counts for missing data.
}

problem.features <- function
### Compute features for one segmentation problem.
(problem.dir
### problemID directory with problemID/coverage.bedGraph.
 ){
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  coverage <- fread(file=file.path(problem.dir, "coverage.bedGraph"))
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

problem.models <- function
### Read models from filesystem. Reads both CSV and RDS files, then
### saves them all in RDS, then deletes all CSV files.
(problem.dir,
### character string path: samples/groupID/sampleID/problems/probID
  verbose=getOption("PeakSegPipeline.verbose", 1)
### print messages?
){
  suffix.vec <- path <- penalty <- penalty.str <- status <-
    NULL
  col.name.list <- list(
    "loss.tsv"=PeakSegDisk::col.name.list$loss,
    "segments.bed"=PeakSegDisk::col.name.list$segments,
    "timing.tsv"=c("penalty", "megabytes", "seconds"))
  colClasses.vec <- c(
    penalty="numeric",
    megabytes="numeric",
    seconds="numeric",
    segments="integer",
    peaks="integer",
    bases="integer",
    bedGraph.lines="integer",
    mean.pen.cost="numeric",
    total.loss="numeric",
    equality.constraints="integer",
    mean.intervals="numeric",
    max.intervals="integer",
    chrom="character",
    chromStart="integer",
    chromEnd="integer",
    status="character",
    mean="numeric")
  csv.file.vec <- Sys.glob(file.path(
    problem.dir, "coverage.bedGraph_penalty=*"))
  csv.err.dt <- if(length(csv.file.vec)==0){
    data.table()
  }else{
    labels.dt <- problem.labels(problem.dir)
    csv.file.dt <- data.table(path=csv.file.vec, nc::capture_first_vec(
      csv.file.vec,
      "penalty=",
      penalty.str=".*?",
      "_",
      suffix.vec=".*"))
    csv.file.dt[, if(
      identical(suffix.vec, names(col.name.list))
    ){
      tsv.data.list <- list()
      for(suffix.i in seq_along(col.name.list)){
        suffix <- suffix.vec[[suffix.i]]
        col.names <- col.name.list[[suffix.i]]
        colClasses <- unname(colClasses.vec[col.names])
        tsv.data.list[[suffix]] <- fread(
          file=path[[suffix.i]],
          col.names=col.names,
          colClasses=colClasses)
      }
      peak.dt <- tsv.data.list$segments.bed[status=="peak"]
      err.dt <- data.table(PeakError::PeakErrorChrom(peak.dt, labels.dt))
      err.row <- with(err.dt, data.table(
        possible.fn=sum(possible.tp),
        possible.fp=sum(possible.fp),
        fn=sum(fn),
        fp=sum(fp)))
      if(nrow(tsv.data.list$timing.tsv)){
        with(tsv.data.list, data.table(
          timing.tsv[loss.tsv, on="penalty"],
          err.row,
          segments.dt=list(list(segments.bed)),
          errors.dt=list(list(err.dt))))
      }
    }, by="penalty.str"][order(penalty)]
  }
  models.rds <- file.path(problem.dir, "models.rds")
  if(file.exists(models.rds)){
    ## If the same penalty.str is found in the rds and csv files, then
    ## take the one from the csv files (more recently computed).
    rds.err.dt <- readRDS(models.rds)
    keep.err.dt <- if(nrow(rds.err.dt)==0){
      data.table()
    }else{
      rds.err.dt[!penalty.str %in% csv.err.dt$penalty.str]
    }
  }else{
    rds.err.dt <- data.table()
    keep.err.dt <- data.table()
  }
  n.rep <- nrow(rds.err.dt)-nrow(keep.err.dt)
  if(verbose)cat(sprintf(
    "Models in csv=%d rds=%d%s\n",
    nrow(csv.err.dt),
    nrow(rds.err.dt),
    if(n.rep)paste0(" replacing ", n.rep) else ""))
  all.err.dt <- rbind(
    if(nrow(csv.err.dt))csv.err.dt else NULL,
    keep.err.dt)
  saveRDS(all.err.dt, models.rds)
  unlink(csv.file.vec)
  all.err.dt
### Data table with all info from model files. One row per model, with
### list columns segments.dt and errors.dt with elements that are data
### tables representing the segmentation model and label errors.
}  

problem.target <- structure(function
### Compute target interval for a segmentation problem. This function
### repeatedly calls PeakSegDisk::PeakSegFPOP_dir with different penalty values,
### until it finds an interval of penalty values with minimal label
### error. The calls to PeakSegFPOP are parallelized using
### future.apply::future_lapply.  A time limit in minutes may be
### specified in a file problem.dir/target.minutes; the search will
### stop at a sub-optimal target interval if this many minutes has
### elapsed. Useful for testing environments with build time limits
### (travis).
(problem.dir,
### problemID directory in which coverage.bedGraph has already been
### computed. If there is a labels.bed file then the number of
### incorrect labels will be computed in order to find the target
### interval of minimal error penalty values.
  verbose=getOption("PeakSegPipeline.verbose", 1)
 ){
  status <- peaks <- errors <- fp <- fn <- penalty <- max.log.lambda <-
    min.log.lambda <- penalty <- . <- done <- total.loss <- mean.pen.cost <-
      bases <- no.next <- is.min <- min.err.interval <- max.lambda <-
        already.computed <- is.other <- dist <- min.lambda <- log.size <-
          mid.lambda <- chrom <- problemStart <- problemEnd <- chromEnd <-
            annotation <- w.fp <- w.fn <- possible.fp <- possible.fn <-
              w.err <- err.min <- is.best <- best.i <- peaks.diff <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  minutes.file <- file.path(problem.dir, "target.minutes")
  minutes.limit <- Inf
  if(file.exists(minutes.file)){
    minutes.dt <- fread(file=minutes.file)
    if(nrow(minutes.dt)==1 && ncol(minutes.dt)==1){
      setnames(minutes.dt, "minutes")
      if(is.numeric(minutes.dt$minutes)){
        minutes.limit <- minutes.dt$minutes
        verbose <- 1 # this happens on travis.
      }
    }
  }
  seconds.start <- as.numeric(Sys.time())
  stopifnot(is.numeric(minutes.limit))
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  problem.coverage(problem.dir)
  labels.dt <- problem.labels(problem.dir)
  if(verbose)cat(nrow(labels.dt), "labels in", problem.dir, "\n")
  model.err.dt <- problem.models(problem.dir)
  error.cols <- c(
    "iteration", "penalty", "fp", "fn", "possible.fp", "possible.fn",
    "peaks", "total.loss")
  ## Compute the label error for one penalty parameter.
  getError <- function(penalty.str){
    stopifnot(is.character(penalty.str))
    stopifnot(length(penalty.str) == 1)
    one.row <- if(penalty.str %in% model.err.dt$penalty.str){
      select.dt <- data.table(penalty.str)
      model.err.dt[select.dt, on="penalty.str"]
    }else{
      result <- PeakSegDisk::PeakSegFPOP_dir(
        problem.dir, penalty.str, problem.tempfile(problem.dir, penalty.str))
      penalty.peaks <- result$segments[status=="peak",]
      tryCatch({
        penalty.error <- PeakErrorChrom(penalty.peaks, labels.dt)
      }, error=function(e)stop(
        "try deleting _segments.bed and recomputing, ",
        "error computing number of incorrect labels: ",
        e)
      )
      with(penalty.error, data.table(
        result$loss,
        possible.fn=sum(possible.tp),
        possible.fp=sum(possible.fp),
        fn=sum(fn),
        fp=sum(fp)))
    }
    data.table(iteration, one.row)[, error.cols, with=FALSE]
  }
  ## Also compute feature vector here so train is faster later.
  problem.features(problem.dir)
  iteration <- 0
  error.list <- if(nrow(model.err.dt)==0){
    list()
  }else{
    split(
      data.table(
        iteration, model.err.dt
      )[, error.cols, with=FALSE],
      model.err.dt$penalty.str)
  }
  last.target.vec <- c(-Inf, Inf)
  target.result.list <- list()
  still.searching <- TRUE
  while(still.searching){
    next.pen <- if(length(error.list)<2){
      c(0, Inf)
    }else{
      error.dt <- do.call(rbind, error.list)[order(penalty)]
      if(!is.numeric(error.dt$penalty)){
        stop("penalty column is not numeric -- check loss in _loss.tsv files")
      }
      error.dt[, errors := fp+fn]
      error.dt[, w.fp := ifelse(possible.fp==0, fp, fp/possible.fp)]
      error.dt[, w.fn := ifelse(possible.fn==0, fn, fn/possible.fn)]
      error.dt[, w.err := w.fp+w.fn]
      unique.peaks <- error.dt[, data.table(
        .SD[which.max(penalty)],
        penalties=.N
      ), by=list(peaks)]
      path.dt <- data.table(penaltyLearning::modelSelection(
        unique.peaks, "total.loss", "peaks"))
      path.dt[, next.pen := max.lambda]
      path.dt[, already.computed := next.pen %in% names(error.list)]
      path.dt[, no.next := c(diff(peaks) == -1, NA)]
      path.dt[, done := already.computed | no.next]
      path.dt[, err.min := errors==min(errors)]
      path.dt[, is.best := FALSE]
      path.dt[err.min==TRUE, is.best := w.err==min(w.err)]
      path.dt[, best.i := cumsum(ifelse(
        c(is.best[1], diff(is.best))==1, 1, 0))]
      if(verbose)print(path.dt[,.(
        penalty, log.pen=log(penalty), peaks, fp, fn, errors, w.err,
        best=ifelse(is.best, best.i, NA))])
      other.candidates <- path.dt[which(0<diff(fn) & diff(fp)<0)]
      interval.dt <- path.dt[is.best==TRUE, {
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
          log.size=max.log.lambda[.N]-min.log.lambda[1],
          peaks.diff=max(peaks)-min(peaks)
        )
      }, by=list(best.i)]
      largest.interval <- interval.dt[which.max(peaks.diff)]
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
      stopping.candidates <- rbind(
        error.candidates, other.candidates)[done==FALSE]
      seconds.now <- as.numeric(Sys.time())
      minutes.elapsed <- (seconds.now-seconds.start)/60
      if(verbose)cat(sprintf(
        "%f minutes elapsed / %f limit\nTarget interval: %f %f change: %f %f\n",
        minutes.elapsed, minutes.limit,
        target.vec[1], target.vec[2],
        diff.target.vec[1], diff.target.vec[2]))
      if(minutes.elapsed < minutes.limit && nrow(stopping.candidates)){
        lambda.vec <- interval.dt[, c(min.lambda, mid.lambda, max.lambda)]
        interval.candidates <- path.dt[next.pen %in% lambda.vec][done==FALSE]
        unique(
          rbind(stopping.candidates, interval.candidates)$next.pen)
      }
    }#if no models else
    if(length(next.pen)==0){
      still.searching <- FALSE
    }else{
      if(verbose)cat(
        "Next =", paste(next.pen, collapse=", "),
        "\n")
      next.str <- paste(next.pen)
      iteration <- iteration+1
      error.list[next.str] <- future.apply::future_lapply(next.str, getError)
    }
  }#while(!is.null(pen))
  list(
    target=target.vec,
    models=problem.models(problem.dir),
    iterations=error.dt)
### List of info related to target interval computation: target is the
### interval of log(penalty) values that achieve minimum incorrect
### labels (numeric vector of length 2), models and iterations are
### data.tables with one row per model.
}, ex=function(){

  library(PeakSegPipeline)
  data(Mono27ac, envir=environment())
  ## Write the Mono27ac data set to disk.
  problem.dir <- file.path(
    tempfile(),
    "H3K27ac-H3K4me3_TDHAM_BP",
    "samples",
    "Mono1_H3K27ac",
    "S001YW_NCMLS",
    "problems",
    "chr11-60000-580000")
  dir.create(problem.dir, recursive=TRUE, showWarnings=FALSE)
  write.table(
    Mono27ac$labels, file.path(problem.dir, "labels.bed"),
    col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
  write.table(
    Mono27ac$coverage, file.path(problem.dir, "coverage.bedGraph"),
    col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

  ## Creating a target.minutes file stops the optimization after that
  ## number of minutes, resulting in an imprecise target interval, but
  ## saving time (to avoid NOTE on CRAN).
  write.table(
    data.frame(minutes=0.05), file.path(problem.dir, "target.minutes"),
    col.names=FALSE, row.names=FALSE, quote=FALSE)

  ## declare future plan for parallel computation.
  if(requireNamespace("future") && interactive()){
    future::plan("multiprocess")
  }

  ## Compute target interval.
  target.list <- problem.target(problem.dir, verbose=1)

  ## These are all the models computed in order to find the target
  ## interval.
  print(target.list$models[order(penalty), list(
    penalty, log.penalty=log(penalty), peaks, total.loss, fn, fp)])

  ## This is the target interval in log(penalty) values.
  print(target.list$target)

})

problem.labels <- function
### read problemID/labels.bed if it exists, otherwise read
### sampleID/labels.bed
(problem.dir
  ## project/samples/groupID/sampleID/problems/problemID
){
  problemStart1 <- problemStart <- chromStart1 <- chromStart <-
    chrom <- problemEnd <- chromEnd <- annotation <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  problem.labels.bed <- file.path(problem.dir, "labels.bed")
  if(file.exists(problem.labels.bed)){
    return(fread(file=problem.labels.bed, col.names=c(
      "chrom", "chromStart", "chromEnd", "annotation")))
  }
  problems.dir <- dirname(problem.dir)
  sample.dir <- dirname(problems.dir)
  sample.labels.bed <- file.path(sample.dir, "labels.bed")
  if(!file.exists(sample.labels.bed)){
    stop(
      "need labels to compute target interval but none found for problem; ",
      "please create either ",
      sample.labels.bed, " or ",
      problem.labels.bed)
  }
  sample.labels <- fread(
    file=sample.labels.bed,
    col.names=c(
      "chrom", "chromStart", "chromEnd", "annotation"))
  if(nrow(sample.labels)==0){
    sample.labels <- data.table(
      chrom=character(),
      chromStart=integer(),
      chromEnd=integer(),
      annotation=character())
  }
  problem <- problem.table(problem.dir)
  problem[, problemStart1 := problemStart +1L]
  sample.labels[, chromStart1 := chromStart +1L]
  setkey(problem, chrom, problemStart1, problemEnd)
  setkey(sample.labels, chrom, chromStart1, chromEnd)
  labels.dt <- foverlaps(problem, sample.labels, nomatch=0L)
  labels.dt[, data.table(
    chrom, chromStart, chromEnd, annotation)]
### data.table with one row for each label and columns chrom,
### chromStart, chromEnd, annotation.
}

problem.predict <- function
### Predict peaks for a genomic segmentation problem.
(problem.dir,
### project/samples/groupID/sampleID/problems/problemID.
  verbose=getOption("PeakSegPipeline.verbose", 1)
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
    warning(
      "Could not compute coverage in", problem.dir,
      "so not predicting peaks.\n")
    return(NULL)
  }
  features.tsv <- file.path(problem.dir, "features.tsv")
  is.computed <- if(file.exists(features.tsv)){
    TRUE
  }else{
    tryCatch({
      problem.features(problem.dir)
      if(verbose)cat(sprintf("Computed %s\n", features.tsv))
      TRUE
    }, error=function(e){
      FALSE
    })
  }
  if(!is.computed){
    if(verbose)cat("Unable to compute", features.tsv, "so not predicting.\n")
    return(NULL)
  }
  features <- fread(file=features.tsv)
  stopifnot(nrow(features)==1)
  feature.mat <- as.matrix(features)
  model.RData <- file.path(data.dir, "model.RData")
  if(file.exists(model.RData)){
    load(model.RData)
  }else{
    stop(
      "Model file ", model.RData,
      " needed for prediction but does not exist; ",
      "run PeakSegPipeline::problem.train('",
      data.dir, "') to create it")
  }
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
  if(verbose)cat(paste0(
    "Predicting penalty=", pred.penalty,
    " log(penalty)=", log(pred.penalty),
    " based on ", n.features,
    " feature", ifelse(n.features==1, "", "s"),
    ".\n"))
  result <- PeakSegDisk::PeakSegFPOP_dir(
    problem.dir, pred.penalty, problem.tempfile(problem.dir, pred.penalty))
  all.peaks <- result$segments[status=="peak", ]
  bases.vec <- all.peaks[, chromEnd-chromStart]
  in.range <- size.model[, lower.bases < bases.vec & bases.vec < upper.bases]
  peaks <- all.peaks[in.range, ]
  ## save peaks.
  peaks.bed <- file.path(problem.dir, "peaks.bed")
  if(verbose)cat(
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
