jobs_create <- function
### Setup a data directory for analysis with PeakSegPipeline.
(data.dir.arg,
### path to project directory.
  verbose=getOption("PeakSegPipeline.verbose", 1)
### TRUE for output, FALSE otherwise.
){
  bases <- problemEnd <- problemStart <- problem.name <- chrom <- row.i <-
    problemStart1 <- chromStart1 <- chromStart <- chromEnd <- chunk.limits <-
      chunk.name <- . <- regions.by.chunk.file <- chunk <- NULL
  ## above to avoid CRAN NOTE.
  if(FALSE){#for debugging.
    data.dir.arg <- "~/genomic-ml/PeakSegFPOP/labels/ATAC_JV_adipose/"
  }
  stop.without.ucsc()
  data.dir <- normalizePath(data.dir.arg, mustWork=TRUE)
  problems.bed <- file.path(data.dir, "problems.bed")
  samples.dir <- file.path(data.dir, "samples")
  problems <- fread(file=problems.bed)
  setnames(problems, c("chrom", "problemStart", "problemEnd"))
  problems[, bases := problemEnd-problemStart]
  problems[, problem.name := sprintf(
    "%s:%d-%d", chrom, problemStart, problemEnd)]
  prob.ord <- problems[order(-bases)]
  prob.ord[, chunk := NA_integer_ ]
  prob.ord[, row.i := 1:.N]
  chunk.id <- 1L
  max.bases.per.chunk <- prob.ord$bases[1]
  big.bases.limit <- max.bases.per.chunk * 0.5
  while(any(is.na(prob.ord$chunk))){
    prob.ord[chunk.id, chunk := chunk.id]
    big.bases <- prob.ord[chunk.id, bases]
    if(big.bases < big.bases.limit){
      small.first <- prob.ord[is.na(chunk)][order(bases)][seq(1, .N, by=30)]
      small.first[, cumsum := big.bases + cumsum(as.numeric(bases))]
      row.vec <- small.first[cumsum<max.bases.per.chunk, row.i]
      prob.ord[row.vec, chunk := chunk.id]
    }
    chunk.id <- chunk.id + 1L
  }
  chunk.stats <- prob.ord[, list(
    problems=.N,
    bases=sum(bases)
  ), by=list(chunk)]
  chunk.stats[, list(
    mean.problems=mean(problems),
    max.problems=max(problems),
    min.bases=min(bases),
    mean.bases=mean(bases),
    max.bases=max(bases)
  )]
  problems[, problemStart1 := problemStart + 1L]
  if(verbose)cat(
    "Read ", nrow(problems),
    " problems from ", problems.bed,
    "\n", sep="")
  unlinkProblems <- function(glob){
    prob.dir.vec <- Sys.glob(glob)
    to.delete <- prob.dir.vec[!basename(prob.dir.vec) %in% problems$problem.name]
    if(length(to.delete)){
      if(verbose)cat("Removing the following old problem directories:\n")
      print(data.table(to.delete))
    }
    unlink(to.delete, recursive=TRUE, force=TRUE)
  }
  sample.dir.glob <- file.path(samples.dir, "*", "*")
  unlinkProblems(file.path(sample.dir.glob, "problems", "*"))
  unlinkProblems(file.path(data.dir, "problems", "*"))
  sample.dir.vec <- Sys.glob(sample.dir.glob)
  labels.bed.vec <- file.path(sample.dir.vec, "labels.bed")
  if(any(file.exists(labels.bed.vec))){
    if(verbose)cat(
      "Some labels.bed files exist, so not running convert_labels\n")
  }else{
    convert_labels(data.dir, verbose=verbose)
  }
  all.job.list <- list()
  for(sample.i in seq_along(sample.dir.vec)){
    sample.dir <- sample.dir.vec[[sample.i]]
    problems.dir <- file.path(sample.dir, "problems")
    labels.bed <- file.path(sample.dir, "labels.bed")
    labels <- if(file.exists(labels.bed))fread(file=labels.bed, col.names=c(
      "chrom", "chromStart", "chromEnd", "annotation"))
    labels.by.problem <- if(length(labels)){
      just.to.check <- PeakError(Peaks(), labels)
      labels[, chromStart1 := chromStart + 1L]
      setkey(labels, chrom, chromStart1, chromEnd)
      setkey(problems, chrom, problemStart1, problemEnd)
      over.dt <- foverlaps(labels, problems, nomatch=0L)
      if(nrow(over.dt) < nrow(labels)){
        warning(
          nrow(labels), " lines in ",
          labels.bed, " but only ",
          nrow(over.dt), " labels occur in ",
          problems.bed)
      }
      split(data.frame(over.dt), over.dt$problem.name)
    }
    if(length(labels.by.problem)){
      if(verbose)cat(sprintf(
        "Step1 %4d / %4d samples %d labels %d labeled problems %s\n",
        sample.i, length(sample.dir.vec),
        nrow(labels), length(labels.by.problem),
        problems.dir))
      all.job.list[[paste("Step1 sample", sample.i)]] <- data.table(
        step=1,
        fun="problem.target",
        arg=file.path(problems.dir, names(labels.by.problem)))
    }
  }#for(sample.i
  all.job.list[["Step2 train"]] <- data.table(
    step=2,
    fun="problem.train",
    arg=data.dir)
  all.job.list[["Step3 predict"]] <- data.table(
    step=3,
    fun="problem.pred.cluster.targets",
    arg=file.path(data.dir, "problems", problems$problem.name))
  all.job.list[["Step4 train joint"]] <- data.table(
    step=4,
    fun="problem.joint.targets.train",
    arg=data.dir)
  all.job.list[["Step5 joint predict"]] <- data.table(
    step=5,
    fun="problem.joint.predict.job",
    arg=file.path(data.dir, "jobs", 1:nrow(problems)))
  all.job.list[["Step6 summarize"]] <- data.table(
    step=6,
    fun="plot_all",
    arg=data.dir)
  all.job <- do.call(rbind, all.job.list)
  hub.sh <- file.path(data.dir, "hub.sh")
  if(!file.exists(hub.sh)){
    code <- sprintf(
      'PeakSegPipeline::create_track_hub("%s", "%s", "%s", "%s")',
      data.dir,
      paste0("http://CHANGE.THIS/~URL/", basename(data.dir)),
      "hg19",
      "toby.hocking@r-project.org")
    script <- sprintf(
      "#!/bin/bash\nRscript -e '%s'",
      code)
    cat(script, file=hub.sh)
  }
  all.job
### data.table with one row for each job and three columns: fun, arg,
### step. fun is the function to call with argument arg, in order
### specified by step (smaller steps first).
}

jobs_submit_batchtools <- structure(function
### Submit PeakSegPipeline jobs via batchtools.
(jobs,
### data.table from jobs_create.
  resources=list(
    walltime = 24*60,#minutes
    memory = 2000,#megabytes per cpu
    ncpus=2,
    ntasks=1,
    chunks.as.arrayjobs=TRUE)
### List of resources for each job, passed to batchtools::submitJobs.
){
  step <- arg <- fun <- NULL
  ## Above to avoid CRAN NOTE.
  requireNamespace("batchtools")
  if(!(
    is.data.table(jobs) &&
    0 < nrow(jobs) &&
    all(c("step", "fun", "arg") %in% names(jobs))
  )){
    stop("jobs must be a data.table with at least one row",
         " and columns step, fun, arg")
  }
  one.job <- jobs[1]
  delete.pattern <- paste0(
    paste(
      "/samples/[^/]+/[^/]+/problems/[^/]+",
      "/problems/[^/]+",
      "/jobs/[^/]+",
      sep="|"),
    "$")
  data.dir <- sub(delete.pattern, "", one.job$arg)
  registry.dir <- file.path(data.dir, "registry")
  dir.create(registry.dir, showWarnings=FALSE)
  steps <- jobs[, list(
    jobs=.N
  ), by=list(step)][order(step)]
  afterok <- NULL
  reg.list <- list()
  for(step.row in 1:nrow(steps)){
    step.i <- steps[step.row, step]
    step.dir <- file.path(registry.dir, step.i)
    unlink(step.dir, recursive=TRUE)
    reg <- reg.list[[paste(step.i)]] <- batchtools::makeRegistry(step.dir)
    step.jobs <- jobs[step==step.i]#[1:min(2, .N)]#for testing
    batchtools::batchMap(function(task.i, job.dt){
      library(PeakSegPipeline)
      job <- job.dt[task.i]
      fun <- get(job$fun)
      fun(job$arg)
    }, 1:nrow(step.jobs), reg=reg, more.args=list(job.dt=step.jobs))
    job.table <- batchtools::getJobTable(reg=reg)
    chunks <- data.table(job.table, chunk=1)
    job.name.vec <- step.jobs[, paste0(
      sub("problem.", "", fun),
      "_",
      basename(arg)
    )]
    batchtools::setJobNames(job.table$job.id, job.name.vec)
    resources[["afterok"]] <- afterok
    batchtools::submitJobs(chunks, resources=resources, reg=reg)
    ##system("squeue -u th798")
    (jobs.done <- batchtools::getJobTable(reg=reg))
    afterok <- sub("_.*", "", jobs.done$batch.id)[[1]]
  }
  reg.list
### A list of registry objects.
}, ex=function(){
  if(FALSE){
    jobs <- jobs_create("~/genomic-ml/PeakSegFPOP/labels/ATAC_JV_adipose/")
    jobs_submit_batchtools(jobs)
  }
})

jobs_create_run <- function
### Run entire PeakSegFPOP + PeakSegJoint pipeline.
(set.dir.path,
### data set directory.
  verbose=getOption("PeakSegPipeline.verbose", 1)
### print messages?
){
  unlink(file.path(
    set.dir.path,
    "samples",
    "*",
    "*",
    "labels.bed"))
  jobs <- jobs_create(set.dir.path, verbose=verbose)
  for(job.i in 1:nrow(jobs)){
    job <- jobs[job.i]
    fun <- get(job$fun)
    fun(job$arg)
  }
}


