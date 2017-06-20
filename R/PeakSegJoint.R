problem.joint.predict.many <- function
### Compute all joint peak predictions for one separate problem.
(prob.dir
### project/problems/problemID
){
  joint.dir.vec <- Sys.glob(file.path(
    prob.dir, "jointProblems", "*"))
  peaks.bed <- file.path(prob.dir, "peaks.bed")
  unlink(peaks.bed)
  prob.progress <- function(joint.dir.i){
    joint.dir <- joint.dir.vec[[joint.dir.i]]
    cat(sprintf(
      "%4d / %4d joint prediction problems %s\n",
      joint.dir.i, length(joint.dir.vec),
      joint.dir))
    jpeaks.bed <- file.path(joint.dir, "peaks.bed")
    already.computed <- if(!file.exists(jpeaks.bed)){
      FALSE
    }else{
      if(0 == file.size(jpeaks.bed)){
        jprob.peaks <- data.table()
        TRUE
      }else{
        tryCatch({
          jprob.peaks <- fread(jpeaks.bed)
          setnames(
            jprob.peaks,
            c("chrom", "chromStart", "chromEnd", "name", "mean"))
          TRUE
        }, error=function(e){
          FALSE
        })
      }
    }
    if(already.computed){
      cat("Skipping since peaks.bed already exists.\n")
    }else{
      jprob.peaks <- problem.joint.predict(joint.dir)
    }
    gc()
    jprob.peaks
  }
  ## out of memory errors, so don't run in parallel!
  peaks.list <- mclapply.or.stop(seq_along(joint.dir.vec), prob.progress)
  ##lapply(seq_along(joint.dir.vec), prob.progress)
  peaks <- if(length(peaks.list)==0){
    data.table()
  }else{
    do.call(rbind, peaks.list)
  }
  ##fread( does not support writing a data.table with 0 rows, so here
  ##we use write.table instead, for convenience.
  write.table(
    peaks,
    peaks.bed,
    quote=FALSE,
    sep="\t",
    col.names=FALSE,
    row.names=FALSE)
  peaks
### data.table of predicted peaks.
}

problem.joint.predict.job <- function
### Compute all joint peak predictions for the joint problems listed
### in jobProblems.bed
(job.dir
### project/jobs/jobID
){
  jobProblems <- fread(file.path(job.dir, "jobProblems.bed"))
  jobs.dir <- dirname(job.dir)
  data.dir <- dirname(jobs.dir)
  problems.dir <- file.path(data.dir, "problems")
  setnames(jobProblems, c("chrom", "problemStart", "problemEnd", "problem.name"))
  jobProblems[, jprob.name := sprintf(
    "%s:%d-%d", chrom, problemStart, problemEnd)]
  prob.progress <- function(joint.dir.i){
    prob <- jobProblems[joint.dir.i]
    joint.dir <- prob[, file.path(
      problems.dir, problem.name, "jointProblems", jprob.name)]
    cat(sprintf(
      "%4d / %4d joint prediction problems %s\n",
      joint.dir.i, nrow(jobProblems),
      joint.dir))
    jpeaks.bed <- file.path(joint.dir, "peaks.bed")
    already.computed <- if(!file.exists(jpeaks.bed)){
      FALSE
    }else{
      if(0 == file.size(jpeaks.bed)){
        jprob.peaks <- data.table()
        loss.dt <- data.table()
        TRUE
      }else{
        tryCatch({
          jprob.peaks <- fread(jpeaks.bed)
          setnames(
            jprob.peaks,
            c("chrom", "chromStart", "chromEnd", "name", "mean"))
          loss.dt <- fread(file.path(joint.dir, "loss.tsv"))
          setnames(loss.dt, "loss.diff")
          TRUE
        }, error=function(e){
          FALSE
        })
      }
    }
    gc()
    jmodel <- if(already.computed){
      cat("Skipping since peaks.bed already exists.\n")
      list(peaks=jprob.peaks, loss=loss.dt)
    }else{
      problem.joint.predict(joint.dir)
    }
    if(nrow(jmodel$peaks)){
      with(jmodel, data.table(
        chrom=peaks$chrom[1],
        peakStart=peaks$chromStart[1],
        peakEnd=peaks$chromEnd[1],
        means=list(peaks[, structure(mean, names=name)]),
        loss.diff=loss$loss.diff,
        problem.name=prob$problem.name
      ))
    }
  }
  jmodel.list <- mclapply.or.stop(1:nrow(jobProblems), prob.progress)
  jobPeaks <- do.call(rbind, jmodel.list)
  jobPeaks.RData <- file.path(job.dir, "jobPeaks.RData")
  save(jobPeaks, file=jobPeaks.RData)
  jobPeaks
### data.table of predicted loss, peak positions, and means per sample
### (in a list column).
}

problem.joint.targets <- function
### Compute targets for a separate problem.
(problem.dir
### project/problems/problemID
 ){
  labels.tsv.vec <- Sys.glob(file.path(
    problem.dir, "jointProblems", "*", "labels.tsv"))
  mclapply.or.stop(seq_along(labels.tsv.vec), function(labels.i){
    labels.tsv <- labels.tsv.vec[[labels.i]]
    jprob.dir <- dirname(labels.tsv)
    cat(sprintf(
      "%4d / %4d labeled joint problems %s\n",
      labels.i, length(labels.tsv.vec),
      jprob.dir))
    target.tsv <- file.path(jprob.dir, "target.tsv")
    if(file.exists(target.tsv)){
      cat("Skipping since target.tsv exists.\n")
    }else{
      problem.joint.target(jprob.dir)
    }
  })
### Nothing.
}

problem.joint.targets.train <- function
### Compute all target intervals then learn a penalty function.
(data.dir
### project directory.
){
  labels.tsv.vec <- Sys.glob(file.path(
    data.dir, "problems", "*", "jointProblems", "*", "labels.tsv"))
  mclapply.or.stop(seq_along(labels.tsv.vec), function(labels.i){
    labels.tsv <- labels.tsv.vec[[labels.i]]
    prob.dir <- dirname(labels.tsv)
    cat(sprintf(
      "%4d / %4d labeled joint problems %s\n",
      labels.i, length(labels.tsv.vec),
      prob.dir))
    target.tsv <- file.path(prob.dir, "target.tsv")
    if(file.exists(target.tsv)){
      cat("Skipping since target.tsv exists.\n")
    }else{
      problem.joint.target(prob.dir)
    }
  })
  problem.joint.train(data.dir)
### Nothing.
}

problem.joint.train <- function
### Learn a penalty function for joint peak prediction.
(data.dir
### project directory.
){
  joint.model.RData <- file.path(data.dir, "joint.model.RData")
  target.tsv.vec <- Sys.glob(file.path(
    data.dir, "problems", "*", "jointProblems", "*", "target.tsv"))
  cat("Found", length(target.tsv.vec), "target.tsv files for training.\n")
  target.mat.list <- list()
  feature.mat.list <- list()
  for(target.tsv.i in seq_along(target.tsv.vec)){
    target.tsv <- target.tsv.vec[[target.tsv.i]]
    target.vec <- scan(target.tsv, quiet=TRUE)
    problem.dir <- dirname(target.tsv)
    segmentations.RData <- file.path(problem.dir, "segmentations.RData")
    load(segmentations.RData)
    if(any(is.finite(target.vec))){
      target.mat.list[[problem.dir]] <- target.vec
      feature.mat.list[[problem.dir]] <- colSums(segmentations$features)
    }
  }
  feature.mat <- do.call(rbind, feature.mat.list)
  target.mat <- do.call(rbind, target.mat.list)
  cat("Training using", nrow(target.mat), "finite targets.\n")
  set.seed(1)
  joint.model <- penaltyLearning::IntervalRegressionCV(
    feature.mat, target.mat,
    min.observations=nrow(feature.mat))
  joint.model$train.mean.vec <- colMeans(feature.mat)
  pred.log.penalty <- joint.model$predict(feature.mat)
  pred.dt <- data.table(
    too.lo=as.logical(pred.log.penalty < target.mat[,1]),
    too.hi=as.logical(target.mat[,2] < pred.log.penalty))
  pred.dt[, status := ifelse(
    too.lo, "low", ifelse(
      too.hi, "high", "correct"))]
  cat("Train errors:\n")
  print(pred.dt[, list(targets=.N), by=status])
  save(joint.model, feature.mat, target.mat, file=joint.model.RData)
  cat("Saved ", joint.model.RData, "\n", sep="")
  jprobs.bed.dt <- data.table(jprobs.bed=Sys.glob(file.path(
    data.dir, "problems", "*", "jointProblems.bed")))
  jprobs <- jprobs.bed.dt[, {
      fread(jprobs.bed)
    }, by=jprobs.bed]
  setnames(jprobs, c(
    "jprobs.bed", "chrom", "problemStart", "problemEnd"))
  jobs.dir <- file.path(data.dir, "jobs")
  job.id.vec <- dir(jobs.dir)
  jprobs[, job := rep(job.id.vec, l=.N) ]
  jprobs[, {
    job.dir <- normalizePath(file.path(jobs.dir, job), mustWork=TRUE)
    fwrite(
      .SD[,.(chrom, problemStart, problemEnd, basename(dirname(jprobs.bed)))],
      file.path(job.dir, "jobProblems.bed"),
      sep="\t", col.names=FALSE)
  }, by=job]
  cat(
    "Saved ",
    length(job.id.vec),
    " jobProblems.bed files to ",
    jobs.dir,
    "\n", sep="")
### Nothing.
}

problem.joint <- function
### Fit a joint model.
(jointProblem.dir
### path/to/jointProblem
){
  segmentations.RData <- file.path(jointProblem.dir, "segmentations.RData")
  problem.bed <- file.path(jointProblem.dir, "problem.bed")
  problem <- fread(problem.bed)
  setnames(problem, c("chrom",  "problemStart", "problemEnd", "problem.name"))
  problem[, problemStart1 := problemStart + 1L]
  setkey(problem, problemStart1, problemEnd)
  jointProblems <- dirname(jointProblem.dir)
  prob.dir <- dirname(jointProblems)
  probs.dir <- dirname(prob.dir)
  data.dir <- dirname(probs.dir)
  samples.dir <- file.path(data.dir, "samples")
  coverage.bedGraph.vec <- Sys.glob(file.path(
    samples.dir, "*", "*", "problems",
    problem$problem.name, "coverage.bedGraph"))
  cat("Found", length(coverage.bedGraph.vec), "samples to jointly segment.\n")
  coverage.list <- list()
  for(coverage.i in seq_along(coverage.bedGraph.vec)){
    coverage.bedGraph <- coverage.bedGraph.vec[[coverage.i]]
    problem.dir <- dirname(coverage.bedGraph)
    save.coverage <- problem[, readCoverage(
      problem.dir, problemStart, problemEnd)]
    coverage.list[[problem.dir]] <- save.coverage
  }
  coverage <- do.call(rbind, coverage.list)
  profile.list <- ProfileList(coverage)
  fit <- PeakSegJointSeveral(coverage)
  segmentations <- ConvertModelList(fit)
  segmentations$features <- featureMatrix(profile.list)
  cat("Writing segmentation and features to", segmentations.RData, "\n")
  save(segmentations, file=segmentations.RData)
  segmentations$coverage <- coverage
  segmentations
### Model from ConvertModelList.
}

readCoverage <- function
### Read sample coverage for one problem from either
### sampleID/coverage.bigWig if it exists, or
### sampleID/problems/problemID/coverage.bedGraph
(problem.dir,
### project/samples/groupID/sampleID/problems/problemID
  start,
### start of coverage to read.
  end
### end of coverage to read.
){
  chrom <- sub(":.*", "", basename(problem.dir))
  jprob <- data.table(chrom, problemStart=start, problemEnd=end)
  problems.dir <- dirname(problem.dir)
  sample.dir <- dirname(problems.dir)
  coverage.bigWig <- file.path(sample.dir, "coverage.bigWig")
  coverage.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
  save.coverage <- if(file.exists(coverage.bigWig)){
    ## If bigWig file has already been computed, then it is much
    ## faster to read it since it is indexed!
    jprob[, readBigWig(
      coverage.bigWig, chrom, problemStart, problemEnd)]
  }else if(0 < file.size(coverage.bedGraph)){#otherwise fread gives error.
    sample.coverage <- fread(
      coverage.bedGraph,
      colClasses=list(NULL=1, integer=2:4))
    setnames(sample.coverage, c("chromStart", "chromEnd", "count"))
    sample.coverage[, chromStart1 := chromStart + 1L]
    setkey(sample.coverage, chromStart1, chromEnd)
    jprob[, problemStart1 := problemStart + 1L]
    setkey(jprob, problemStart1, problemEnd)
    problem.coverage <- foverlaps(sample.coverage, jprob, nomatch=0L)
    problem.coverage[chromStart < problemStart, chromStart := problemStart]
    problem.coverage[problemEnd < chromEnd, chromEnd := problemEnd]
    problem.coverage[chromStart < chromEnd,]
    ## Note that problem.coverage (from coverage.bedGraph) is
    ## guaranteed to have rows with 0 coverage -- this is required for
    ## PeakSegFPOP! but readBigWig really just returns what is stored
    ## in the bigWig -- if there are no rows with 0 coverage, then
    ## there will be none in the output save.coverage. This is no
    ## problem for running the PeakSegJoint model, but it could be a
    ## problem for plotting (error in geom_step if only one row to
    ## plot).
  }
  ## If we don't have the if() below, we get Error in
  ## data.table(sample.id, sample.group,
  ## problem.coverage[chromStart < : Item 3 has no length.
  if(is.data.table(save.coverage) && 0 < nrow(save.coverage)){
    sample.id <- basename(sample.dir)
    group.dir <- dirname(sample.dir)
    sample.group <- basename(group.dir)
    data.table(sample.id, sample.group, save.coverage)
  }
### Either the data.table of coverage, or NULL if no coverage data exists.
}

problem.joint.predict <- function
### Compute peak predictions for a joint problem.
(jointProblem.dir
### project/problems/problemID/jointProblems/jointProbID
){
  segmentations.RData <- file.path(jointProblem.dir, "segmentations.RData")
  if(file.exists(segmentations.RData)){
    cat("Loading model from ", segmentations.RData, "\n", sep="")
    load(segmentations.RData)
  }else{
    segmentations <- problem.joint(jointProblem.dir)
  }
  chrom <- sub(":.*", "", basename(jointProblem.dir))
  jprobs.dir <- dirname(jointProblem.dir)
  prob.dir <- dirname(jprobs.dir)
  probs.dir <- dirname(prob.dir)
  set.dir <- dirname(probs.dir)
  joint.model.RData <- file.path(set.dir, "joint.model.RData")
  load(joint.model.RData)
  feature.mat <- rbind(colSums(segmentations$features))
  stopifnot(nrow(feature.mat)==1)
  if(length(feature.mat)==length(joint.model$train.mean.vec)){
    is.bad <- !is.finite(feature.mat)
    feature.mat[is.bad] <- joint.model$train.mean.vec[is.bad]
  }else{
    stop(
      "feature.mat has ", length(feature.mat),
      " columns but ", length(joint.model$train.mean.vec),
      " features were used to train the joint model")
  }    
  log.penalty <- as.numeric(joint.model$predict(feature.mat))
  stopifnot(length(log.penalty)==1)
  stopifnot(is.finite(log.penalty))
  selected <- subset(
    segmentations$modelSelection,
    min.log.lambda < log.penalty & log.penalty < max.log.lambda)
  loss.tsv <- file.path(jointProblem.dir, "loss.tsv")
  if(selected$peaks == 0){
    unlink(loss.tsv)
    pred.dt <- data.table()
    loss.dt <- data.table()
  }else{
    selected.loss <- segmentations$loss[paste(selected$peaks), "loss"]
    flat.loss <- segmentations$loss["0", "loss"]
    loss.dt <- data.table(
      loss.diff=flat.loss-selected.loss)
    write.table(
      loss.dt,
      loss.tsv,
      quote=FALSE,
      sep="\t",
      col.names=FALSE,
      row.names=FALSE)
    pred.df <- subset(segmentations$peaks, peaks==selected$peaks)
    pred.dt <- with(pred.df, data.table(
      chrom,
      chromStart,
      chromEnd,
      name=paste0(sample.group, "/", sample.id),
      mean))
  }
  peaks.bed <- file.path(jointProblem.dir, "peaks.bed")
  cat("Writing ",
      nrow(pred.dt), " peaks to ",
      peaks.bed,
      "\n", sep="")
  write.table(
    pred.dt, peaks.bed,
    quote=FALSE,
    sep="\t",
    col.names=FALSE,
    row.names=FALSE)
  list(peaks=pred.dt, loss=loss.dt)
### list of peaks and loss.
}

problem.joint.target <- function
### Compute target interval for a joint problem.
(jointProblem.dir
### Joint problem directory.
){
  segmentations.RData <- file.path(jointProblem.dir, "segmentations.RData")
  if(file.exists(segmentations.RData)){
    cat("Loading model from ", segmentations.RData, "\n", sep="")
    load(segmentations.RData)
  }else{
    segmentations <- problem.joint(jointProblem.dir)
  }
  labels.bed <- file.path(jointProblem.dir, "labels.tsv")
  labels <- fread(labels.bed)
  setnames(labels, c(
    "chrom", "chromStart", "chromEnd", "annotation",
    "sample.id", "sample.group"))
  jprob <- fread(file.path(jointProblem.dir, "problem.bed"))
  setnames(jprob, c("chrom", "problemStart", "problemEnd", "problem.name"))
  ##   [peakStart]   [peakEnd]   labels
  ## ______   ___________  ______ joint problems
  ## ____________________________
  ## peakStart must end inside the joint problem,
  ## and peakEnd must start inside the joint problem.
  ## otherwise the annotation should be considered noPeaks.
  labels[{
    (annotation=="peakEnd" & chromStart < jprob$problemStart) |
      (annotation=="peakStart" & jprob$problemEnd < chromEnd)
  }, annotation := "noPeaks"]
  fit.error <- PeakSegJointError(segmentations, labels)
  if(FALSE){
    show.peaks <- 8
    show.peaks.df <- subset(segmentations$peaks, peaks==show.peaks)
    show.errors <- fit.error$error.regions[[paste(show.peaks)]]
    ann.colors <-
      c(noPeaks="#f6f4bf",
        peakStart="#ffafaf",
        peakEnd="#ff4c4c",
        peaks="#a445ee")
    ggplot()+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "lines"))+
      facet_grid(sample.group + sample.id ~ ., scales="free")+
      scale_fill_manual(values=ann.colors)+
      geom_tallrect(aes(
        xmin=chromStart/1e3,
        xmax=chromEnd/1e3,
        fill=annotation),
        color="grey",
        alpha=0.5,
        data=labels)+
      geom_tallrect(aes(
        xmin=chromStart/1e3,
        xmax=chromEnd/1e3,
        linetype=status),
        color="black",
        size=1,
        fill=NA,
        data=show.errors)+
      scale_linetype_manual(
        "error type",
        limits=c("correct", 
                 "false negative",
                 "false positive"),
        values=c(correct=0,
                 "false negative"=3,
                 "false positive"=1))+
      geom_step(aes(chromStart/1e3, count),
                data=segmentations$coverage,
                color="grey50")+
      geom_segment(aes(chromStart/1e3, 0,
                       xend=chromEnd/1e3, yend=0),
                   data=show.peaks.df,
                   color="deepskyblue",
                   size=2)
    ## geom_segment(aes(chromStart/1e3, mean,
    ##                  xend=chromEnd/1e3, yend=mean),
    ##              data=
    ##              color="green")
  }
  cat("Train error:\n")
  print(fit.error$modelSelection[, c(
    "min.log.lambda", "max.log.lambda", "peaks", "errors")])
  target.tsv <- file.path(jointProblem.dir, "target.tsv")
  cat(
    "Writing target interval (",
    paste(fit.error$target, collapse=", "),
    ") to ", 
    target.tsv,
    "\n", sep="")
  write(fit.error$target, target.tsv, sep="\t")
### list of output from PeakSegJointError.
}

problem.joint.plot <- function
### Plot one chunk.
(chunk.dir
### project/problems/problemID/chunks/chunkID
){
  if(!require(ggplot2)){
    stop("please install ggplot2")
  }
  chunks.dir <- dirname(chunk.dir)
  prob.dir <- dirname(chunks.dir)
  prob.name <- basename(prob.dir)
  probs.dir <- dirname(prob.dir)
  proj.dir <- dirname(probs.dir)
  problem.dir.vec <- Sys.glob(file.path(
    proj.dir, "samples", "*", "*", "problems", prob.name))
  chunk <- fread(file.path(chunk.dir, "chunk.bed"))
  setnames(chunk, c("chrom", "chunkStart", "chunkEnd"))
  chunk[, chunkStart1 := chunkStart + 1L]
  setkey(chunk, chunkStart1, chunkEnd)
  labels <- fread(file.path(chunk.dir, "labels.tsv"))
  cat("Read",
      nrow(labels),
      "labels.\n")
  jointProblems <- fread(file.path(prob.dir, "jointProblems.bed"))
  setnames(jointProblems, c("chrom", "problemStart", "problemEnd"))
  jointProblems[, problemStart1 := problemStart + 1L]
  jointProblems[, problem.name := sprintf(
    "%s:%d-%d", chrom, problemStart, problemEnd)]
  setkey(jointProblems, problemStart1, problemEnd)
  probs.in.chunk <- foverlaps(jointProblems, chunk, nomatch=0L)
  probs.in.chunk$sample.group <- "problems"
  probs.in.chunk$sample.id <- "joint"
  cat("Read",
      nrow(jointProblems),
      "joint problems, plotting",
      nrow(probs.in.chunk),
      "in chunk.\n")
  coverage.list <- list()
  separate.peaks.list <- list()
  for(sample.i in seq_along(problem.dir.vec)){
    problem.dir <- problem.dir.vec[[sample.i]]
    problems.dir <- dirname(problem.dir)
    sample.dir <- dirname(problems.dir)
    sample.id <- basename(sample.dir)
    group.dir <- dirname(sample.dir)
    sample.group <- basename(group.dir)
    chunk.cov <- chunk[, readCoverage(problem.dir, chunkStart, chunkEnd)]
    coverage.list[[problem.dir]] <- chunk.cov
    ## Also store peaks in this chunk, if there are any.
    sample.peaks <- tryCatch({
      fread(file.path(problem.dir, "peaks.bed"))
    }, error=function(e){
      data.table()
    })
    if(nrow(sample.peaks)){
      setnames(sample.peaks, c("chrom", "peakStart", "peakEnd", "status", "mean"))
      sample.peaks[, peakStart1 := peakStart + 1L]
      setkey(sample.peaks, peakStart1, peakEnd)
      chunk.peaks <- foverlaps(sample.peaks, chunk, nomatch=0L)
      if(nrow(chunk.peaks)){
        separate.peaks.list[[problem.dir]] <- data.table(
          sample.id, sample.group, chunk.peaks)
      }
    }
  }
  coverage <- do.call(rbind, coverage.list)
  cat("Read",
      length(coverage.list),
      "samples of coverage.\n")
  cat("Read",
      length(separate.peaks.list),
      "samples of separate peak predictions.\n")
  joint.peaks.list <- list()
  for(joint.i in 1:nrow(probs.in.chunk)){
    prob <- probs.in.chunk[joint.i,]
    tryCatch({
      peaks <- fread(file.path(
        prob.dir, "jointProblems", prob$problem.name, "peaks.bed"))
      setnames(peaks, c("chrom", "peakStart", "peakEnd", "sample.path", "mean"))
      peaks[, sample.id := sub(".*/", "", sample.path)]
      peaks[, sample.group := sub("/.*", "", sample.path)]
      joint.peaks.list[[prob$problem.name]] <- peaks
    }, error=function(e){
      NULL
    })
  }
  cat("Read",
      length(joint.peaks.list),
      "joint peak predictions.\n")
  joint.peaks <- do.call(rbind, joint.peaks.list)
  ann.colors <-
    c(noPeaks="#f6f4bf",
      peakStart="#ffafaf",
      peakEnd="#ff4c4c",
      peaks="#a445ee")
  gg <- ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(sample.group + sample.id ~ ., scales="free")+
    scale_y_continuous(
      "aligned read coverage",
      breaks=function(limits){
        lim <- floor(limits[2])
        if(lim==0){
          Inf
        }else{
          lim
        }
      })+
    scale_x_continuous(paste(
      "position on",
      coverage$chrom[1],
      "(kb = kilo bases)"))+
    ## geom_tallrect(aes(
    ##   xmin=problemStart/1e3,
    ##   xmax=problemEnd/1e3),
    ##   alpha=0.5,
    ##   size=3,
    ##   color="black",
    ##   fill=NA,
    ##   data=probs.in.chunk)+
    geom_segment(aes(
      problemStart/1e3, 0,
      xend=problemEnd/1e3, yend=0),
      size=1,
      color="blue",
      data=probs.in.chunk)+
    geom_point(aes(
      problemStart/1e3, 0),
      color="blue",
      data=probs.in.chunk)+
    geom_tallrect(aes(
      xmin=chromStart/1e3, 
      xmax=chromEnd/1e3,
      fill=annotation), 
      alpha=0.5,
      data=labels)+
    scale_fill_manual("label", values=ann.colors)+
    scale_color_manual(values=c(separate="black", joint="deepskyblue"))+
    scale_size_manual(values=c(separate=2, joint=3))+
    geom_rect(aes(
      xmin=chromStart/1e3, xmax=chromEnd/1e3,
      ymin=0, ymax=count),
      data=coverage,
      color="grey50")
  if(length(joint.peaks)){
    joint.peaks$peak.type <- "joint"
    gg <- gg+
      geom_point(aes(
        peakStart/1e3, 0,
        color=peak.type,
        size=peak.type),
        data=joint.peaks)+
      geom_segment(aes(
        peakStart/1e3, 0,
        xend=peakEnd/1e3, yend=0,
        color=peak.type,
        size=peak.type),
        data=joint.peaks)
  }
  if(length(separate.peaks.list)){
    separate.peaks <- do.call(rbind, separate.peaks.list)
    separate.peaks$peak.type <- "separate"
    gg <- gg+
      geom_segment(aes(
        peakStart/1e3, 0,
        xend=peakEnd/1e3, yend=0,
        color=peak.type,
        size=peak.type),
                   data=separate.peaks)+
      geom_point(aes(
        peakStart/1e3, 0,
        color=peak.type,
        size=peak.type),
                 data=separate.peaks)
  }
  n.rows <- length(coverage.list) + 2
  mypng <- function(base, g){
    f <- file.path(chunk.dir, base)
    cat("Writing ",
        f,
        "\n", sep="")
    cairo.limit <- 32767
    h <- 60*n.rows
    if(cairo.limit < h){
      h <- floor(cairo.limit / n.rows) * n.rows
    }
    png(f, res=100, width=1000, height=h)
    print(g)
    dev.off()
    thumb.png <- sub(".png$", "-thumb.png", f)
    cmd <- sprintf("convert %s -resize 230 %s", f, thumb.png)
    system(cmd)
  }
  mypng("figure-predictions-zoomout.png", gg)
  gg.zoom <- gg+
    coord_cartesian(
      xlim=chunk[, c(chunkStart, chunkEnd)/1e3],
      expand=FALSE)
  mypng("figure-predictions.png", gg.zoom)
### Nothing
}  
