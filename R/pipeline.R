pipeline <- function
### Run entire PeakSegFPOP + PeakSegJoint pipeline.
(set.dir.path,
### data set directory.
  verbose=0
### print messages?
){
  set.dir <- normalizePath(set.dir.path, mustWork=TRUE)
  ## First convert labels.
  convert_labels(set.dir)
  ## Create problems for each sample.
  create_problems_all(set.dir)
  ## Compute target interval for each problem.
  samples.dir <- file.path(set.dir, "samples")
  labels.bed.vec <- Sys.glob(file.path(
    samples.dir, "*", "*", "problems", "*", "labels.bed"))
  lapply(labels.bed.vec, function(labels.bed){
    sample.dir <- dirname(labels.bed)
    problem.target(sample.dir, verbose=verbose)
  })
  ## Train single-sample model.
  problem.train(set.dir)
  ## Single-sample prediction and peak clustering, one job for each
  ## problem.
  problem.dir.vec <- Sys.glob(file.path(set.dir, "problems", "*"))
  for(problem.dir in problem.dir.vec){
    problem.name <- basename(problem.dir)
    create.cmd <- paste(
      "bash",
      file.path(problem.dir, "jointProblems.bed.sh"))
    ## This includes lapply in parallel on samples.
    system.or.stop(create.cmd)
  }
  ## Compute target intervals for multi-sample problems, then learn a
  ## penalty function for joint peak prediction.
  problem.joint.targets.train(set.dir)
  ## Joint prediction, new job scripts.
  job.sh.vec <- Sys.glob(file.path(
    set.dir, "jobs", "*", "jobPeaks.sh"))
  for(job.sh in job.sh.vec){
    job.cmd <- paste("bash", job.sh)
    system.or.stop(job.cmd)
  }
  ## Summarize peak predictions on a web page.
  peaks.tsv.sh <- file.path(set.dir, "peaks_matrix.tsv.sh")
  final.cmd <- paste("bash", peaks.tsv.sh)
  system.or.stop(final.cmd)
}
