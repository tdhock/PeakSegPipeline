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
  jprob.name <- chrom <- problemStart <- problemEnd <- problem.name <-
    jprob.name <- sample.loss.diff <- group.loss.diff <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
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
    peakInfo.rds <- file.path(joint.dir, "peakInfo.rds")
    already.computed <- if(!file.exists(peakInfo.rds)){
      FALSE
    }else{
      tryCatch({
        pred.row <- readRDS(peakInfo.rds)
        TRUE
      }, error=function(e){
        FALSE
      })
    }
    if(already.computed){
      cat("Skipping since peakInfo.rds already exists.\n")
    }else{
      pred.row <- problem.joint.predict(joint.dir)
    }
    if(pred.row[, sample.loss.diff==0 && group.loss.diff==0]){
      data.table()
    }else{
      prob[, data.table(problem.name, jprob.name, pred.row)]
    }
  }
  jmodel.list <- mclapply.or.stop(1:nrow(jobProblems), prob.progress)
  jobPeaks <- do.call(rbind, jmodel.list)
  jobPeaks.RData <- file.path(job.dir, "jobPeaks.RData")
  save(jobPeaks, file=jobPeaks.RData)
  jobPeaks
### data.table of predicted peaks, one row for each job, same columns
### as from problem.joint.predict.
}

problem.joint.targets <- function
### Compute targets for a separate problem.
(problem.dir
### project/problems/problemID
 ){
  segmentations <- model <- min.log.lambda <- max.log.lambda <- NULL
  ## above variable defined in RData file.
  labels.tsv.vec <- Sys.glob(file.path(
    problem.dir, "jointProblems", "*", "labels.tsv"))
  targets.features.list <- mclapply.or.stop(
    seq_along(labels.tsv.vec), function(labels.i){
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
      target.dt <- fread(target.tsv)
      setkey(target.dt, model)
      segmentations.RData <- file.path(jprob.dir, "segmentations.RData")
      load(segmentations.RData)
      data.table(
        sample.target=list(
          target.dt["sample", c(min.log.lambda, max.log.lambda)]),
        group.target=list(
          target.dt["group", c(min.log.lambda, max.log.lambda)]),
        features=list(colSums(segmentations$features)))
    })
  targets.features <- do.call(rbind, targets.features.list)
  jointTargets.list <- lapply(targets.features, function(L)do.call(rbind, L))
  jointTargets.rds <- file.path(problem.dir, "jointTargets.rds")
  saveRDS(jointTargets.list, jointTargets.rds)
  jointTargets.list
### Named list of two matrices: targets and features.
}

problem.joint.targets.train <- function
### Compute all joint target intervals then learn joint penalty
### functions.
(data.dir
### project directory.
){
  problem.dir.vec <- Sys.glob(file.path(
    data.dir, "problems", "*"))
  mclapply.or.stop(seq_along(problem.dir.vec), function(problem.i){
    problem.dir <- problem.dir.vec[[problem.i]]
    cat(sprintf(
      "%4d / %4d problems %s\n",
      problem.i, length(problem.dir.vec),
      problem.dir))
    jointTargets.rds <- file.path(problem.dir, "jointTargets.rds")
    if(file.exists(jointTargets.rds)){
      cat("Skipping since jointTargets.rds exists.\n")
    }else{
      problem.joint.targets(problem.dir)
    }
  })
  problem.joint.train(data.dir)
### Nothing.
}

problem.pred.cluster.targets <- function
### For a given problem, predict independently for each sample, then
### cluster peaks across samples, then compute joint target intervals.
(proj.dir
### project directory.
){
  peaks.dt <- problem.predict.allSamples(proj.dir)
  create_problems_joint(proj.dir, peaks.dt)
  problem.joint.targets(proj.dir)
### List of features and target matrices (same as problem.joint.target).
}

problem.joint.train <- function
### Learn a penalty function for joint peak prediction.
(data.dir
### project directory.
){
  segmentations <- status <- too.lo <- too.hi <- jprobs.bed <- job <-
    . <- chrom <- problemStart <- problemEnd <- jointTargets.rds <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  joint.model.RData <- file.path(data.dir, "joint.model.RData")
  jointTargets.rds.vec <- Sys.glob(file.path(
    data.dir, "problems", "*", "jointTargets.rds"))
  ## Below we read the features and targets into a data.table with
  ## list columns.
  jointDT <- data.table(jointTargets.rds=jointTargets.rds.vec)[, {
      L <- readRDS(jointTargets.rds)
      lapply(L, list)
  }, by=list(jointTargets.rds)]
  mat.list <- lapply(jointDT[,-1,with=FALSE], function(L)do.call(rbind, L))
  model.list <- list()
  for(model.name in c("sample", "group")){
    target.name <- paste0(model.name, ".target")
    all.target.mat <- mat.list[[target.name]]
    keep <- apply(is.finite(all.target.mat), 1, any)
    n.keep <- sum(keep)
    if(n.keep<2){
      stop("need at least two targets to train model")
    }
    cat("Training", model.name, "model using", n.keep, "finite targets.\n")
    target.mat <- all.target.mat[keep,]
    feature.mat <- mat.list$features[keep,]
    set.seed(1)
    joint.model <- penaltyLearning::IntervalRegressionCV(
      feature.mat, target.mat,
      min.observations=n.keep)
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
    model.list[[model.name]] <- joint.model
  }
  save(model.list, mat.list, file=joint.model.RData)
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
  problemStart1 <- problemStart <- problemEnd <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
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
  if(FALSE){
    ggplot()+
      geom_rect(aes(
        xmin=chromStart/1e3, xmax=chromEnd/1e3,
        ymin=0, ymax=count),
        data=coverage)+
      facet_grid(sample.id + sample.group ~ .)+
      theme_bw()+
      theme(panel.spacing=grid::unit(0, "lines"))
  }
  profile.list <- PeakSegJoint::ProfileList(coverage)
  segmentations <- PeakSegJoint::PeakSegJointFaster(profile.list)
  segmentations$features <- PeakSegJoint::featureMatrixJoint(profile.list)
  cat("Writing segmentation and features to", segmentations.RData, "\n")
  save(segmentations, file=segmentations.RData)
  segmentations$coverage <- coverage
  segmentations
### Model from PeakSegJointFaster.
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
  problemStart <- problemEnd <- chromStart1 <- chromStart <-
    chromEnd <- problemStart1 <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  chrom <- sub(":.*", "", basename(problem.dir))
  jprob <- data.table(chrom, problemStart=start, problemEnd=end)
  problems.dir <- dirname(problem.dir)
  sample.dir <- dirname(problems.dir)
  coverage.bigWig <- file.path(sample.dir, "coverage.bigWig")
  save.coverage <- jprob[, readBigWig(
    coverage.bigWig, chrom, problemStart, problemEnd)]
  ## If we don't have the if() below, we get Error in
  ## data.table(sample.id, sample.group,
  ## problem.coverage[chromStart < : Item 3 has no length.
  if(nrow(save.coverage)==0){
    save.coverage <- data.table(chromStart=start, chromEnd=end, count=0L)
  }
  sample.id <- basename(sample.dir)
  group.path <- dirname(sample.dir)
  sample.group <- basename(group.path)
  data.table(sample.id, sample.group, save.coverage)
### The data.table of coverage, which can be 1 row with 0 coverage if
### no coverage data exists. (this is necessary for standard output of
### the right number of samples for each joint problem, even those
### with no coverage in some samples)
}

problem.joint.predict <- function
### Compute peak predictions for a joint problem.
(jointProblem.dir
### project/problems/problemID/jointProblems/jointProbID
){
  joint.model <- min.log.lambda <- max.log.lambda <- peaks <- NULL
  model.list <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  segmentations.RData <- file.path(jointProblem.dir, "segmentations.RData")
  if(file.exists(segmentations.RData)){
    cat("Loading model from ", segmentations.RData, "\n", sep="")
    load(segmentations.RData)
  }else{
    segmentations <- problem.joint(jointProblem.dir)
  }
  joint.prob <- basename(jointProblem.dir)
  chrom <- sub(":.*", "", joint.prob)
  jprobs.dir <- dirname(jointProblem.dir)
  prob.dir <- dirname(jprobs.dir)
  probs.dir <- dirname(prob.dir)
  set.dir <- dirname(probs.dir)
  joint.model.RData <- file.path(set.dir, "joint.model.RData")
  objs <- load(joint.model.RData)
  feature.mat <- rbind(colSums(segmentations$features))
  stopifnot(nrow(feature.mat)==1)
  bkg.vec <- with(segmentations, rowMeans(mean_mat[,c(1,3)]))
  bkg.vec[!segmentations$is.feasible] <- NA
  peak.str <- with(segmentations, {
    sprintf("%s:%d-%d", chrom, peak_start_end[1], peak_start_end[2])
  })
  listmat <- function(x, N=names(x)){
    list(Matrix(x, length(N), 1, dimnames=list(
      N, peak.str)))
  }
  pred.dt <- with(segmentations, data.table(
    ## Need to save mean values of all peaks,
    ## to compute peak height matrix later.
    chrom,
    peakStart=peak_start_end[1],
    peakEnd=peak_start_end[2],
    peak.mean.vec=listmat(mean_mat[,2]),
    sample.loss.diff.vec=listmat(flat_loss_vec-peak_loss_vec, sample.id),
    ## Also save mean values of background which is lower than the peak mean,
    ## to compute peak height matrix.
    background.mean.vec=listmat(bkg.vec),
    ## sample.peaks and group.peaks start out all FALSE
    ## but below we will set some entries to TRUE
    ## based on where we predict peaks.
    sample.peaks.vec=listmat(FALSE, sample.id),
    group.peaks.vec=listmat(FALSE, names(group.list))))
  for(model.name in names(model.list)){
    joint.model <- model.list[[model.name]]
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
    m <- function(suffix)paste0(model.name, suffix)
    selected <- subset(
      segmentations[[m(".modelSelection")]],
      min.log.lambda < log.penalty & log.penalty < max.log.lambda)
    loss.diff <- if(selected$complexity == 0){
      0
    }else{
      with.peaks <- segmentations[[m(".loss.diff.vec")]][1:selected$complexity]
      pred.dt[[m(".peaks.vec")]][[1]][names(with.peaks), ] <- TRUE
      -sum(with.peaks)
    }
    pred.dt[[paste0(model.name, ".loss.diff")]] <- loss.diff
  }
  peakInfo.rds <- file.path(jointProblem.dir, "peakInfo.rds")
  cat("Writing ", peakInfo.rds, "\n", sep="")
  saveRDS(pred.dt, peakInfo.rds)
  pred.dt
### data.table with one row, describing predicted peaks for both group
### and sample models. Columns are chrom, peakStart, peakEnd,
### peak.mean.vec, background.mean.vec (for computing normalized peak
### height relative to background), sample.loss.diff.vec (for ranking
### peaks in a given sample from most to least likely), sample.peaks,
### group.peaks, sample.loss.diff, group.loss.diff (for ranking
### likelihood of peak regions in joint model).
}

problem.joint.target <- function
### Compute target interval for a joint problem.
(jointProblem.dir
### Joint problem directory.
){
  annotation <- chromStart <- chromEnd <- sample.path <-
    sample.group <- sample.id <- flat.errors <- peak.errors <-
      errors <- complexity <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
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
  labels[, sample.path := paste0(sample.group, "/", sample.id)]
  peaks.dt <- data.table(
    sample.path=names(segmentations$sample.loss.diff.vec))
  setkey(peaks.dt, sample.path)
  ## we assign chromStart/end after to avoid error "Item 3 has no length"
  peaks.dt[, chromStart := segmentations$peak_start_end[1] ]
  peaks.dt[, chromEnd := segmentations$peak_start_end[2] ]
  group.dt <- with(segmentations, {
    data.table(sample.group=names(group.list))[, list(
      sample.path=group.list[[sample.group]]
    ), by=sample.group]
  })
  setkey(group.dt, sample.path)
  peak.errors.dt <- labels[, {
    sample.peaks <- peaks.dt[sample.path, nomatch=0L]
    error.df <- PeakError::PeakErrorChrom(sample.peaks, .SD)
    with(error.df, data.table(peak.errors=sum(fp+fn)))
  }, by=sample.path]
  setkey(peak.errors.dt, sample.path)
  flat.errors.dt <- labels[, {
    error.df <- PeakError::PeakErrorChrom(PeakError::Peaks(), .SD)
    with(error.df, data.table(flat.errors=sum(fp+fn)))
  }, by=sample.path]
  setkey(flat.errors.dt, sample.path)
  errors.all.samples <- flat.errors.dt[peak.errors.dt[group.dt]]
  errors.all.samples[is.na(flat.errors), flat.errors := 0]
  errors.all.samples[is.na(peak.errors), peak.errors := 0]
  errors.selected.samples <-
    errors.all.samples[names(segmentations$sample.loss.diff.vec)]
  errors.all.groups <- errors.all.samples[, list(
    flat.errors=sum(flat.errors),
    peak.errors=sum(peak.errors)
  ), by=sample.group]
  errors.selected.groups <-
    errors.all.groups[names(segmentations$group.loss.diff.vec)]
  flat.errors.total <- sum(flat.errors.dt$flat.errors)
  sample.errors.vec <- flat.errors.total + errors.selected.samples[, c(
    0, cumsum(peak.errors)-cumsum(flat.errors))]
  sample.select <- data.table(segmentations$sample.modelSelection)
  sample.select[, errors := sample.errors.vec[complexity+1] ]
  group.errors.vec <- flat.errors.total + errors.selected.groups[, c(
    0, cumsum(peak.errors)-cumsum(flat.errors))]
  group.select <- data.table(segmentations$group.modelSelection)
  group.select[, errors := group.errors.vec[complexity+1] ]
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
      theme(panel.spacing=grid::unit(0, "lines"))+
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
  cat("Train error for samples:\n")
  print(sample.select[, c(
    "min.log.lambda", "max.log.lambda", "complexity", "errors")])
  cat("Train error for groups:\n")
  print(group.select[, c(
    "min.log.lambda", "max.log.lambda", "complexity", "errors")])
  both.select <- rbind(
    data.table(sample.select, model="sample"),
    data.table(group.select, model="group"))
  target.dt <- penaltyLearning::targetIntervals(both.select, "model")
  cat("Target intervals of minimum error penalty values:\n")
  print(target.dt)
  target.tsv <- file.path(jointProblem.dir, "target.tsv")
  cat(
    "Writing target intervals to ", 
    target.tsv,
    "\n", sep="")
  fwrite(target.dt, target.tsv)
### Nothing.
}

problem.joint.plot <- function
### Plot one chunk.
(chunk.dir
### project/problems/problemID/chunks/chunkID
){
  chunkStart1 <- chunkStart <- chunkEnd <- problemStart1 <-
    problem.name <- chrom <- problemEnd <- problemStart <- peakStart1 <-
      peakStart <- peakEnd <- sample.path <- chromStart <- chromEnd <-
        annotation <- count <- peak.type <- NULL
  sample.peaks.vec <- peak.mean.vec <- NULL
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
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
    peakInfo <- readRDS(file.path(
      prob.dir, "jointProblems", prob$problem.name, "peakInfo.rds"))
    joint.peaks.list[[prob$problem.name]] <- peakInfo[, {
      is.peak <- as.logical(sample.peaks.vec[[1]])
      if(any(is.peak)){
        mean.vec <- peak.mean.vec[[1]]
        sample.path <- rownames(mean.vec)[is.peak]
        data.table(
          chrom, peakStart, peakEnd,
          sample.path,
          mean=mean.vec[is.peak],
          sample.id=sub(".*/", "", sample.path),
          sample.group=sub("/.*", "", sample.path))
      }
    }]
  }
  joint.peaks <- do.call(rbind, joint.peaks.list)
  cat(
    "Read joint peak predictions:",
    nrow(joint.peaks), "peaks in",
    length(joint.peaks.list), "genomic regions,",
    nrow(probs.in.chunk), "peakInfo.RData files.\n"
  )
  ann.colors <-
    c(noPeaks="#f6f4bf",
      peakStart="#ffafaf",
      peakEnd="#ff4c4c",
      peaks="#a445ee")
  gg <- ggplot()+
    theme_bw()+
    theme(panel.spacing=grid::unit(0, "lines"))+
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
    scale_fill_manual(
      "label", values=ann.colors)+
    scale_color_manual(
      values=c(separate="black", joint="deepskyblue"))+
    scale_size_manual(
      values=c(separate=2, joint=3))+
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
  ##mypng("figure-predictions-zoomout.png", gg)
  gg.zoom <- gg+
    coord_cartesian(
      xlim=chunk[, c(chunkStart, chunkEnd)/1e3],
      expand=FALSE)
  mypng("figure-predictions.png", gg.zoom)
### Nothing
}  
