### Setup a data directory for analysis with PeakSegPipeline.
jobs_create <- function(data.dir.arg, verbose=FALSE){
  if(FALSE){
    data.dir.arg <- "~/genomic-ml/PeakSegFPOP/labels/ATAC_JV_adipose/"
  }
  data.dir <- normalizePath(data.dir.arg, mustWork=TRUE)
  problems.bed <- file.path(data.dir, "problems.bed")
  samples.dir <- file.path(data.dir, "samples")
  problems <- fread(problems.bed)
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
  all.job.list <- list()
  for(sample.i in seq_along(sample.dir.vec)){
    sample.dir <- sample.dir.vec[[sample.i]]
    problems.dir <- file.path(sample.dir, "problems")
    labels.bed <- file.path(sample.dir, "labels.bed")
    labels <- tryCatch({
      labels <- fread(labels.bed)
      setnames(labels, c("chrom", "chromStart", "chromEnd", "annotation"))
      labels
    }, error=function(e){
      if(verbose)cat("No labels in", labels.bed, "\n")
      data.table()
    })
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
    if(verbose)cat(sprintf(
      "Writing %4d / %4d samples %d labels %d labeled problems %s\n",
      sample.i, length(sample.dir.vec),
      nrow(labels), length(labels.by.problem),
      problems.dir))
    all.job.list[[paste("Step1 sample", sample.i)]] <- data.table(
      step=1,
      fun="problem.target",
      arg=file.path(problems.dir, names(labels.by.problem)))
    ## all.job.list[[paste("Step3 predict", sample.i)]] <- data.table(
    ##   step=3,
    ##   fun="problem.predict",
    ##   arg=file.path(problems.dir, problems$problem.name))
  }#for(sample.i
  chunk.limits.RData <- file.path(data.dir, "chunk.limits.RData")
  if(file.exists(chunk.limits.RData)){
    objs <- load(chunk.limits.RData)
    chunks <- data.table(chunk.limits)
    chunks[, chunk.name := sprintf("%s:%d-%d", chrom, chromStart, chromEnd)]
    chunks[, chromStart1 := chromStart+1L]
    setkey(chunks, chrom, chromStart1, chromEnd)
    chunks.with.problems <- foverlaps(problems, chunks, nomatch=0L)
    setkey(chunks.with.problems, problem.name)
  }
  ## Now write data_dir/problems/*/jointProblems.bed.sh
  all.job.list[["Step2 train"]] <- data.table(
    step=2,
    fun="problem.train",
    arg=data.dir)
  all.job.list[["Step3 predict"]] <- data.table(
    step=3,
    fun="problem.pred.cluster.targets",
    arg=file.path(data.dir, "problems", problems$problem.name))
  for(problem.i in 1:nrow(problems)){
    problem <- problems[problem.i,]
    prob.dir <- file.path(data.dir, "problems", problem$problem.name)
    dir.create(prob.dir, showWarnings=FALSE, recursive=TRUE)
    if(file.exists(chunk.limits.RData) &&
       problem$problem.name %in% chunks.with.problems$problem.name){
      ## write a directory for every chunk.
      problem.chunks <- chunks.with.problems[problem$problem.name]
      if(verbose)cat("Writing ", nrow(problem.chunks),
          " chunks in ", prob.dir,
          "\n", sep="")
      for(chunk.i in seq_along(problem.chunks$chunk.name)){
        chunk <- problem.chunks[chunk.i,]
        chunk.dir <- file.path(prob.dir, "chunks", chunk$chunk.name)
        dir.create(chunk.dir, showWarnings=FALSE, recursive=TRUE)
        chunk.bed <- file.path(chunk.dir, "chunk.bed")
        fwrite(
          chunk[, .(chrom, chromStart, chromEnd)],
          chunk.bed,
          sep="\t",
          col.names=FALSE,
          quote=FALSE)
        chunk.labels <- regions.by.chunk.file[[paste(chunk$file.and.chunk)]]
        labels.tsv <- file.path(chunk.dir, "labels.tsv")
        fwrite(
          chunk.labels,
          labels.tsv,
          sep="\t",
          col.names=TRUE,
          quote=FALSE)
      }
    }
    ## joint prediction jobs script.
    job.dir <- file.path(data.dir, "jobs", problem.i)
    dir.create(job.dir, showWarnings=FALSE, recursive=TRUE)
  }#for(problem.i
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
  if(FALSE){
    all.job[, basename := basename(arg)]
    all.job[, problem.name := ifelse(
      grepl(":", basename),
      basename,
      NA_character_)]
    prob.chunk <- prob.ord[, list(chunk, problem.name)]
    job.chunks <- prob.chunk[all.job, on=list(problem.name)]
    job.chunks[is.na(chunk), chunk := 1:.N, by=list(step)]
    job.chunks[, list(jobs=.N, chunks=length(unique(chunk))), by=list(step)]
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
    walltime = 3600, memory = 1024, ncpus=1, ntasks=1,
    chunks.as.arrayjobs=TRUE)
### List of resources for each job, passed to batchtools::submitJobs.
){
  requireNamespace("batchtools")
  data.dir <- jobs[step==2, arg]
  registry.dir <- file.path(data.dir, "registry")
  dir.create(registry.dir, showWarnings=FALSE)
  steps <- jobs[, list(
    jobs=.N
  ), by=list(step)][order(step)]
  afterok <- NULL
  for(step.i in 1:nrow(steps)){
    step.dir <- file.path(registry.dir, step.i)
    unlink(step.dir, recursive=TRUE)
    reg <- batchtools::makeRegistry(step.dir)
    step.jobs <- jobs[step==step.i]#[1:min(2, .N)]#for testing
    batchtools::batchMap(function(task.i, job.dt){
      library(PeakSegPipeline)
      job <- job.dt[task.i]
      fun <- get(job$fun)
      fun(job$arg)
    }, 1:nrow(step.jobs), reg=reg, more.args=list(job.dt=step.jobs))
    job.table <- batchtools::getJobTable(reg=reg)
    chunks <- data.table(job.table, chunk=1)
    batchtools::setJobNames(job.table$job.id, step.jobs[, paste(
      fun, basename(arg))])
    resources[["afterok"]] <- afterok
    batchtools::submitJobs(chunks, resources=resources, reg=reg)
    ##system("squeue -u th798")
    (jobs.done <- batchtools::getJobTable(reg=reg))
    afterok <- sub("_.*", "", jobs.done$batch.id)[[1]]
  }
}, ex=function(){
  if(FALSE){
    jobs <- jobs_create("~/genomic-ml/PeakSegFPOP/labels/ATAC_JV_adipose/")
    jobs_submit_batchtools(jobs)
  }
})


