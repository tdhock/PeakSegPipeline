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

PeakSegFPOP_disk <- structure(function
### Run the PeakSeg Functional Pruning Optimal Partitioning algorithm,
### using a file on disk (rather than in memory as in
### coseg::PeakSegFPOP) to store the O(N) function piece lists, each of size O(log N).
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
  names.list <- PeakSegFPOP_disk(tmp, "10.5")
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
  peaks.list <- mclapply.or.stop(problem.vec, problem.predict)
  do.call(rbind, peaks.list)
### data.table of predicted peaks.
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
  chrom <- problemStart <- problemEnd <- count.num.str <- coverage <-
    count.int <- count.int.str <- chromStart <- chromEnd <- J <-
      count <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
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
    timing <- fread(penalty_timing.tsv)
    setnames(timing, c("penalty", "megabytes", "seconds"))
    first.line <- fread(paste("head -1", penalty_segments.bed))
    setnames(first.line, c("chrom", "chromStart", "chromEnd", "status", "mean"))
    last.line <- fread(paste("tail -1", penalty_segments.bed))
    setnames(last.line, c("chrom", "chromStart", "chromEnd", "status", "mean"))
    penalty.loss <- fread(penalty_loss.tsv)
    setnames(penalty.loss, c(
      "penalty", "segments", "peaks", "bases",
      "mean.pen.cost", "total.cost", "status",
      "mean.intervals", "max.intervals"))
    loss.segments.consistent <-
      first.line$chromEnd-last.line$chromStart == penalty.loss$bases
    pattern <- paste0(
      "(?<chrom>chr[^:]+)",
      ":",
      "(?<chromStart>[0-9]+)",
      "-",
      "(?<chromEnd>[0-9]+)")
    problem <- str_match_named(basename(problem.dir), pattern, list(
      chromStart=as.integer,
      chromEnd=as.integer))
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
### Compute target interval for a segmentation problem. This function
### repeated calls problem.PeakSegFPOP with different penalty values,
### until it finds an interval of penalty values with minimal label
### error. The calls to PeakSegFPOP are parallelized using mclapply if
### you set options(mc.cores).
(problem.dir
### problemID directory in which coverage.bedGraph has already been
### computed. If there is a labels.bed file then the number of
### incorrect labels will be computed in order to find the target
### interval of minimal error penalty values.
 ){
  status <- peaks <- errors <- fp <- fn <- penalty <- max.log.lambda <-
    min.log.lambda <- penalty <- . <- done <- total.cost <- mean.pen.cost <-
      bases <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
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
  ## Compute the target interval given the errors computed in dt.
  getTarget <- function(dt){
    peaks.tab <- table(dt$peaks)
    error.sorted <- dt[order(peaks), ][c(TRUE, diff(peaks) != 0),]
    error.sorted[, errors := fp + fn]
    setkey(error.sorted, peaks)
    ##error.sorted[, model.complexity := oracleModelComplexity(bases, segments)]
    path <- penaltyLearning::modelSelection(
      error.sorted, "total.cost", "peaks")
    path.dt <- data.table(path)
    setkey(path.dt, peaks)
    join.dt <- error.sorted[path.dt][order(penalty),]
    direction.list <- list(start=1, end=-1)
    side.vec.list <- list(fn="end", fp="start", errors=c("start", "end"))
    result <- list(models=path, candidates=list())
    for(error.col in c("fp", "fn", "errors")){
      indices <- penaltyLearning::largestContinuousMinimumC(
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
  iteration <- 0
  target.result.list <- list()
  while(length(next.pen)){
    cat(
      "Next =", paste(next.pen, collapse=", "),
      "mc.cores=", getOption("mc.cores"),
      "\n")
    next.str <- paste(next.pen)
    iteration <- iteration+1
    error.list[next.str] <- mclapply.or.stop(next.str, getError)
    error.dt <- do.call(rbind, error.list)[order(-penalty),]
    if(!is.numeric(error.dt$penalty)){
      stop("penalty column is not numeric -- check loss in _loss.tsv files")
    }
    print(error.dt[,.(penalty, peaks, status, fp, fn)])
    target.list <- getTarget(error.dt)
    target.vec <- c(
      target.list$candidates[["errors start"]]$min.log.lambda,
      target.list$candidates[["errors end"]]$max.log.lambda)
    target.result.list[[paste(iteration)]] <- data.table(
      iteration,
      min.log.lambda=target.vec[1],
      max.log.lambda=target.vec[2])
    is.error <- grepl("error", names(target.list$candidates))
    error.candidates <- do.call(rbind, target.list$candidates[is.error])
    other.candidates <- do.call(rbind, target.list$candidates[!is.error])
    other.in.target <- other.candidates[done==FALSE &
        target.vec[1] < log(next.pen) & log(next.pen) < target.vec[2],]
    next.pen <- if(nrow(other.in.target)==0){
      cat("No fp/fn min in min(error) interval => refine limits of min(error) interval.\n")
      ##print(error.candidates)
      error.candidates[done==FALSE, unique(next.pen)]
    }else{      
      ## Inside the minimum error interval, we have found a spot where
      ## the fn or fp reaches a minimum. This means that we should try
      ## exploring a few penalty values between the fp/fn limits.
      pen.vec <- other.candidates[done==FALSE, sort(unique(next.pen))]
      ##print(other.candidates)
      print(pen.vec)
      if(length(pen.vec)==1){
        ## There is only one unique value, so explore it. This is
        ## possible for an error profile of 3 2 1 1 2 3............
        ## (fp = 3 2 1 0 0 0, fn = 0 0 0 1 2 3) in which case we just
        ## want to explore between the ones.
        cat("FP/FN min in min(error) interval and one unique value to explore.\n")
        pen.vec
      }else{
        ## Rather than simply evaluating the penalties at the borders,
        ## we try a grid of penalties on the log scale.
        cat("FP/FN min in min(error) interval and several values to explore.\n")
        exp(seq(log(pen.vec[1]), log(pen.vec[2]), l=4))
      }
    }
    if(FALSE && interactive() && length(next.pen) && requireNamespace("ggplot2Animint")){
      ## This may cause a crash if executed within mclapply, and it is
      ## mostly just for debugging purposes.
      gg <- ggplot2Animint::ggplot()+
        ggplot2Animint::geom_abline(
          ggplot2Animint::aes(slope=peaks, intercept=total.cost),
                    data=error.dt)+
        ggplot2Animint::geom_vline(
          ggplot2Animint::aes(xintercept=penalty),
                   color="red",
                   data=data.table(penalty=next.pen))+
        ggplot2Animint::geom_point(
          ggplot2Animint::aes(penalty, mean.pen.cost*bases),
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
    target.iterations=do.call(rbind, target.result.list),
    models=error.dt)
### List of info related to target interval computation: target is the
### interval of log(penalty) values that achieve minimum incorrect
### labels (numeric vector of length 2), target.iterations is a
### data.table with target intervals as a function of iteration,
### models is a data.table with one row per model for which the label
### error was computed.
}

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
