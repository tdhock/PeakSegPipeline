problem.train <- function(data.dir.str){
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
    features.list[[problem.dir]] <- fread(features.tsv)
    targets.list[[problem.dir]] <- scan(target.tsv, quiet=TRUE)
  }
  features <- as.matrix(do.call(rbind, features.list))
  targets <- do.call(rbind, targets.list)
  set.seed(1)
  model <- if(nrow(features) < 10){
    IntervalRegressionUnregularized(
      features[, c("log.quartile.100%", "log.data")], targets)
  }else{
    IntervalRegressionCV(
      features, targets, verbose=0,
      initial.regularization=1e-4,
      min.observations=nrow(features),
      reg.type=ifelse(nrow(features) < 20, "1sd", "min(mean)"))
  }
  model$train.mean.vec <- colMeans(features)
  cat("Learned regularization parameter and weights:\n")
  print(model$pred.param.mat)
  pred.log.penalty <- as.numeric(model$predict(features))
  pred.dt <- data.table(
    problem.dir=dirname(target.tsv.vec),
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
  pred.dt[, list(targets=.N), by=status]
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
  size.plot <- ggplot()+
    geom_histogram(
      aes(
        log10.bases, ..density..),
      data=correct.peaks)+
    geom_vline(
      aes(xintercept=mean),
      data=size.model,
      size=1, color="red")+
    PeakSegJoint::geom_tallrect(
      aes(xmin=lower.lim, xmax=upper.lim),
      data=size.model, fill="red")+
    geom_line(
      aes(log10.bases, prob),
      data=normal.dens, color="red", size=1)+
    geom_text(
      aes(log10.bases, 0, label=label, hjust=hjust),
      data=base.labels, vjust=1)
  cat("Writing model to", model.RData, "\n")
  save(
    model, features, targets,
    size.model, size.plot, correct.peaks,
    file=model.RData)
}

PeakSegFPOP_disk <- structure(function
### Run the PeakSeg Functional Pruning Optimal Partitioning algorithm,
### using a file on disk (rather than in memory as in
### coseg::PeakSegFPOP) to store the function piece lists.
### Finds the optimal change-points using the Poisson loss and the
### PeakSeg constraint. For N data points, the functional pruning
### algorithm is O(N log N) time and disk space, and O(log N) memory.
### It computes the exact
### solution to the following optimization problem. Let Z be an
### N-vector of count data (count.vec, non-negative integers), let W
### be an N-vector of positive weights (weight.vec), and let penalty
### be a non-negative real number. Find the N-vector M of real numbers
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
### columns: chrom, chromStart, chromEnd, coverage.
  pen.str
### character scalar that can be converted to a numeric scalar via
### as.numeric: non-negative penalty. More penalty means fewer
### peaks. 0 and Inf are OK.
){
  if(!(
    is.character(bedGraph.file) &&
    length(bedGraph.file)==1 &&
    file.exists(bedGraph.file)
  )){
    stop("bedGraph.file must be the name of a data file to segment")
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
    PACKAGE="PeakSegPipeline")
  prefix <- paste0(bedGraph.file, "_penalty=", pen.str)
  result$segments <- paste0(prefix, "_segments.bed")
  result$db <- paste0(prefix, ".db")
  result$loss <- paste0(prefix, "_loss.tsv")
  result
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
  names.list <- PeakSegFPOP_disk(tmp, 10)
  unlink(names.list$db)
  seg.df <- read.table(names.list$segments)
  names(seg.df) <- c("chrom", "chromStart", "chromEnd", "status", "mean")
  seg.df
  loss.df <- read.table(names.list$loss)
  names(loss.df) <- c(
    "penalty", "segments", "peaks", "bases",
    "mean.pen.cost", "total.cost", "status",
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
  mclapply.or.stop(problem.vec, problem.predict)
### List of data tables (predicted peaks).
}

problem.coverage <- function
### Ensure that coverage.bedGraph has been correctly computed for a
### particular genomic segmentation problem.
(problem.dir
### Path to a directory like sampleID/problems/problemID where
### sampleID/coverage.bedGraph contains counts of aligned reads in the
### entire genome, and sampleID/problems/problemID/problem.bed
### contains one line that indicates the genomic coordinates of a
### particular segmentation problem.
 ){
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  problems.dir <- dirname(problem.dir)
  sample.dir <- dirname(problems.dir)
  problem.bed <- file.path(problem.dir, "problem.bed")
  problem <- fread(problem.bed)
  setnames(problem, c("chrom", "problemStart", "problemEnd"))
  ## First check if problem/coverage.bedGraph has been created.
  prob.cov.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
  coverage.ok <- tryCatch({
    tail.cmd <- paste("tail -1", prob.cov.bedGraph)
    last.cov <- fread(tail.cmd)
    setnames(last.cov, c("chrom", "chromStart", "chromEnd", "coverage"))
    last.cov$chromEnd == problem$problemEnd
  }, error=function(e){
    FALSE
  })
  ## If problemID/coverage.bedGraph has already been computed, than we
  ## have nothing to do.
  if(!coverage.ok){
    ## Create problemID/coverage.bedGraph from either
    ## sampleID/coverage.bigWig or sampleID/coverage.bedGraph.
    coverage.bedGraph <- file.path(sample.dir, "coverage.bedGraph")
    coverage.bigWig <- file.path(sample.dir, "coverage.bigWig")
    cov.cmd <- if(file.exists(coverage.bigWig)){
      problem[, sprintf(
        "bigWigToBedGraph -chrom=%s -start=%d -end=%d %s %s",
        chrom, problemStart, problemEnd,
        coverage.bigWig, prob.cov.bedGraph)]
    }else if(file.exists(coverage.bedGraph)){
      ## Use intersectBed -sorted to avoid memory
      ## problems. coverage.bedGraph needs to be -a since that is
      ## reported in the output.
      paste(
        "intersectBed -sorted",
        "-a", coverage.bedGraph,
        "-b", problem.bed,
        ">", prob.cov.bedGraph)
    }else{
      stop("To compute ", prob.cov.bedGraph,
           " need either ", coverage.bigWig,
           " or ", coverage.bedGraph,
           " which do not exist.")
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
    min.above.zero <- prob.cov[0 < coverage, min(coverage)]
    prob.cov[, count.num := coverage/min.above.zero]
    prob.cov[, count.num.str := paste(count.num)]
    prob.cov[, count.int := as.integer(round(count.num))]
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
    fwrite(
      zero.cov,
      prob.cov.bedGraph,
      quote=FALSE,
      sep="\t",
      col.names=FALSE)
  }
### Nothing. If necessary, the intersectBed or bigWigToBedGraph
### command line program is used to create problemID/coverage.bedGraph
### and then we (1) stop if there are any negative data, (2) stop if
### the data are not integers, or can not be normalized by the
### smallest non-zero value to obtain integers, and (3) add lines with
### zero counts for missing data.
}

problem.PeakSegFPOP <- function
### Run PeakSegFPOP on one genomic segmentation problem directory.
(problem.dir,
### Path to a directory like sampleID/problems/problemID which
### contains a coverage.bedGraph file with the aligned read counts for
### one genomic segmentation problem.
 penalty.str
### Penalty parameter to pass to the PeakSegFPOP command line program.
 ){
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  stopifnot(is.character(penalty.str))
  stopifnot(length(penalty.str)==1)
  prob.cov.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
  pre <- paste0(prob.cov.bedGraph, "_penalty=", penalty.str)
  penalty_segments.bed <- paste0(pre, "_segments.bed")
  penalty_loss.tsv <- paste0(pre, "_loss.tsv")
  penalty_timing.tsv <- paste0(pre, "_timing.tsv")
  already.computed <- tryCatch({
    first.line <- fread(paste("head -1", penalty_segments.bed))
    setnames(first.line, c("chrom", "chromStart", "chromEnd", "status", "mean"))
    last.line <- fread(paste("tail -1", penalty_segments.bed))
    setnames(last.line, c("chrom", "chromStart", "chromEnd", "status", "mean"))
    penalty.loss <- fread(penalty_loss.tsv)
    setnames(penalty.loss, c(
      "penalty", "segments", "peaks", "bases",
      "mean.pen.cost", "total.cost", "status",
      "mean.intervals", "max.intervals"))
    if(first.line$chromEnd-last.line$chromStart == penalty.loss$bases){
      TRUE
    }else{
      FALSE
    }
  }, error=function(e){
    FALSE
  })
  if(already.computed){
    timing <- fread(penalty_timing.tsv)
    setnames(timing, c("penalty", "megabytes", "seconds"))
  }else{
    penalty.db <- paste0(pre, ".db")
    seconds <- system.time({
      PeakSegFPOP_disk(prob.cov.bedGraph, penalty.str)
    })[["elapsed"]]
    megabytes <- if(file.exists(penalty.db)){
      file.size(penalty.db)/1024/1024
    }else{
      0
    }
    timing <- data.table(
      penalty.str,
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
      "mean.pen.cost", "total.cost", "status",
      "mean.intervals", "max.intervals"))
  }
  penalty.segs <- fread(penalty_segments.bed)
  setnames(penalty.segs, c("chrom","chromStart", "chromEnd", "status", "mean"))
  list(
    segments=penalty.segs,
    loss=penalty.loss,
    timing=timing)
### List of data.tables: segments has one row for every segment in the
### optimal model, loss has one row and contains the Poisson loss and
### feasibility, and timing is one row with the time and disk usage.
}  

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

problem.target <- function
### Compute target interval for a segmentation problem.
(problem.dir
### problemID directory in which coverage.bedGraph has already been
### computed. If there is a labels.bed file then the number of
### incorrect labels will be computed in order to find the target
### interval of minimal error penalty values.
 ){
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  c.info <- problem.coverage(problem.dir)
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
    penalty.error <- PeakErrorChrom(penalty.peaks, problem.labels)
    with(penalty.error, data.table(
      result$loss,
      fn=sum(fn),
      fp=sum(fp)))
  }  
  ## Compute the target interval given the errors computed in dt.
  getTarget <- function(dt){
    peaks.tab <- table(dt$peaks)
    error.sorted <- dt[order(peaks), ][c(TRUE, diff(peaks) != 0),]
    error.sorted[, errors := fp + fn]
    setkey(error.sorted, peaks)
    ##error.sorted[, model.complexity := oracleModelComplexity(bases, segments)]
    path <- error.sorted[, penaltyLearning::modelSelection(
      total.cost, peaks, peaks)]
    path.dt <- data.table(path)
    setkey(path.dt, peaks)
    join.dt <- error.sorted[path.dt][order(penalty),]
    direction.list <- list(start=1, end=-1)
    side.vec.list <- list(fn="end", fp="start", errors=c("start", "end"))
    result <- list(models=path, candidates=list())
    for(error.col in c("fp", "fn", "errors")){
      indices <- largestContinuousMinimum(
        join.dt[[error.col]],
        join.dt[, max.log.lambda-min.log.lambda]
        )
      side.vec <- side.vec.list[[error.col]]
      for(side in side.vec){
        direction <- direction.list[[side]]
        index <- indices[[side]]
        model <- join.dt[index,]
        index.outside <- index - direction
        neighboring.peaks <- model$peaks + direction
        found.neighbor <- neighboring.peaks %in% join.dt$peaks
        multiple.penalties <- if(index.outside %in% seq_along(join.dt$peaks)){
          model.outside <- join.dt[index.outside,]
          peaks.num <- c(model.outside$peaks, model$peaks)
          peaks.str <- paste(peaks.num)
          peaks.counts <- peaks.tab[peaks.str]
          any(1 < peaks.counts)
        }else{
          FALSE
        }
        ## cost + lambda * model.complexity =
        ## cost + penalty * peaks =>
        ## penalty = lambda * model.complexity / peaks.
        ## lambda is output by exactModelSelection,
        ## penalty is input by PeakSegFPOP.
        next.pen <- ifelse(side=="start", model$min.lambda, model$max.lambda)
        already.computed <- paste(next.pen) %in% names(error.list)
        done <- found.neighbor | multiple.penalties | already.computed
        result$candidates[[paste(error.col, side)]] <- data.table(
          model, found.neighbor, multiple.penalties, already.computed,
          done, next.pen)
      }
    }
    result
  }
  error.list <- list()
  next.pen <- c(0, Inf)
  while(length(next.pen)){
    cat("Next =", paste(next.pen, collapse=", "), "\n")
    next.str <- paste(next.pen)
    error.list[next.str] <- mclapply.or.stop(next.str, getError)
    error.dt <- do.call(rbind, error.list)[order(-penalty),]
    print(error.dt[,.(penalty, peaks, status, fp, fn)])
    target.list <- getTarget(error.dt)
    target.vec <- c(
      target.list$candidates[["errors start"]]$min.log.lambda,
      target.list$candidates[["errors end"]]$max.log.lambda)
    is.error <- grepl("error", names(target.list$candidates))
    error.candidates <- do.call(rbind, target.list$candidates[is.error])
    other.candidates <- do.call(rbind, target.list$candidates[!is.error])
    other.in.target <- other.candidates[done==FALSE &
        target.vec[1] < log(next.pen) & log(next.pen) < target.vec[2],]
    next.pen <- if(nrow(other.in.target)){
      ## Inside the minimum error interval, we have found a spot where
      ## the fn or fp reaches a minimum. This means that we should try
      ## exploring a few penalty values between the fp/fn limits.
      pen.vec <- other.candidates[done==FALSE, sort(unique(next.pen))]
      if(length(pen.vec)==1){
        ## There is only one unique value, so explore it. This is
        ## possible for an error profile of 3 2 1 1 2 3............
        ## (fp = 3 2 1 0 0 0, fn = 0 0 0 1 2 3) in which case we just
        ## want to explore between the ones.
        pen.vec
      }else{
        ## Rather than simply evaluating the penalties at the borders,
        ## we try a grid of penalties on the log scale.
        exp(seq(log(pen.vec[1]), log(pen.vec[2]), l=4))
      }
    }else{
      error.candidates[done==FALSE, unique(next.pen)]
    }
    if(interactive() && length(next.pen)){
      gg <- ggplot()+
        geom_abline(aes(slope=peaks, intercept=total.cost),
                    data=error.dt)+
        geom_vline(aes(xintercept=penalty),
                   color="red",
                   data=data.table(penalty=next.pen))+
        geom_point(aes(penalty, mean.pen.cost*bases),
                   data=error.dt)
      print(gg)
    }
  }#while(!is.null(pen))
  write.table(
    error.dt,
    file.path(problem.dir, "target_models.tsv"),
    sep="\t",
    quote=FALSE,
    row.names=FALSE,
    col.names=TRUE)
  write(target.vec, file.path(problem.dir, "target.tsv"), sep="\t")
  ## Also compute feature vector here so train is faster later.
  problem.features(problem.dir)
  list(
    target=target.vec,
    models=error.dt)
### List of info related to target interval computation: target is the
### interval of log(penalty) values that achieve minimum incorrect
### labels (numeric vector of length 2), models is a data.table with
### one row per model for which the label error was computed, selected
### is a data.frame that describes which penalty values will select
### which models.
}

problem.predict <- function
### Predict peaks for a genomic segmentation problem.
(problem.dir
### project/samples/groupID/sampleID/problems/problemID.
 ){
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  problems.dir <- dirname(problem.dir)
  sample.dir <- dirname(problems.dir)
  group.dir <- dirname(sample.dir)
  samples.dir <- dirname(group.dir)
  data.dir <- dirname(samples.dir)
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
  peaks
### data.table of peak predictions.
}
