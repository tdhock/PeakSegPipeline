## Edit the following definition to reflect your cluster
## configuration.
PBS.header <- "#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=24:00:00
#PBS -A bws-221-ae
#PBS -m ae
#PBS -M tdhock5@gmail.com
#PBS -V"

## input: a directory with labels.bed and coverage.bedGraph.

## output: problems directory with sub-directories for each problem
## with labels.

arg.vec <- "test/input"
arg.vec <- commandArgs(trailingOnly=TRUE)
if(length(arg.vec) != 1){
  stop("usage: Rscript create_problems_all.R data_dir")
}

library(data.table)
library(PeakError)

Rscript <- function(...){
  code <- sprintf(...)
  if(any(grepl("'", code))){
    print(code)
    stop("there can not be any ' in code")
  }
  sprintf("Rscript -e '%s'", code)
}

data.dir <- normalizePath(arg.vec[1], mustWork=TRUE)
problems.bed <- file.path(data.dir, "problems.bed")
samples.dir <- file.path(data.dir, "samples")

problems <- fread(problems.bed)
setnames(problems, c("chrom", "problemStart", "problemEnd"))
problems[, problemStart1 := problemStart + 1L]
problems[, problem.name := sprintf(
  "%s:%d-%d", chrom, problemStart, problemEnd)]
cat(
  "Read ", nrow(problems),
  " problems from ", problems.bed,
  "\n", sep="")
unlinkProblems <- function(glob){
  prob.dir.vec <- Sys.glob(glob)
  to.delete <- prob.dir.vec[!basename(prob.dir.vec) %in% problems$problem.name]
  if(length(to.delete)){
    cat("Removing the following old problem directories:\n")
    print(data.table(to.delete))
  }
  unlink(to.delete, recursive=TRUE, force=TRUE)
}
sample.dir.glob <- file.path(samples.dir, "*", "*")
unlinkProblems(file.path(sample.dir.glob, "problems", "*"))
unlinkProblems(file.path(data.dir, "problems", "*"))

sample.dir.vec <- Sys.glob(sample.dir.glob)
for(sample.i in seq_along(sample.dir.vec)){
  sample.dir <- sample.dir.vec[[sample.i]]
  problems.dir <- file.path(sample.dir, "problems")
  coverage.bedGraph <- file.path(sample.dir, "coverage.bedGraph")
  labels.bed <- file.path(sample.dir, "labels.bed")
  labels <- tryCatch({
    labels <- fread(labels.bed)
    setnames(labels, c("chrom", "chromStart", "chromEnd", "annotation"))
    labels
  }, error=function(e){
    cat("No labels in", labels.bed, "\n")
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
  makeProblem <- function(problem.i){
    problem <- data.frame(problems)[problem.i,]
    problem.dir <- file.path(problems.dir, problem$problem.name)
    ## cat(sprintf(
    ##   "%4d / %4d %s\n",
    ##   problem.i, nrow(problems), problem.dir))
    dir.create(problem.dir, showWarnings=FALSE, recursive=TRUE)
    problem.bed <- file.path(problem.dir, "problem.bed")
    prob.text <- with(problem, data.frame(
      chrom,
      chromStart=sprintf("%d", problemStart),
      chromEnd=sprintf("%d", problemEnd)))
    write.table(
      prob.text, problem.bed,
      quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    if(problem$problem.name %in% names(labels.by.problem)){
      problem.labels <- labels.by.problem[[problem$problem.name]]
      prob.lab.bed <- file.path(problem.dir, "labels.bed")
      lab.text <- with(problem.labels, data.frame(
        chrom,
        chromStart=sprintf("%d", chromStart),
        chromEnd=sprintf("%d", chromEnd),
        annotation))
      write.table(
        lab.text,
        prob.lab.bed,
        quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    }
    ## Script for coverage.
    target.tsv <- file.path(problem.dir, "target.tsv")
    sh.file <- paste0(target.tsv, ".sh")
    target.cmd <- Rscript('coseg::problem.target("%s")', problem.dir)
    script.txt <- paste0(PBS.header, "
#PBS -o ", target.tsv, ".out
#PBS -e ", target.tsv, ".err
#PBS -N Target", problem$problem.name, "
", target.cmd, "
")
    writeLines(script.txt, sh.file)
    ## Script for peaks.
    peaks.bed <- file.path(problem.dir, "peaks.bed")
    sh.file <- paste0(peaks.bed, ".sh")
    predict.cmd <- Rscript(
      'coseg::problem.predict("%s")',
      problem.dir)
    script.txt <- paste0(PBS.header, "
#PBS -o ", peaks.bed, ".out
#PBS -e ", peaks.bed, ".err
#PBS -N Predict", problem$problem.name, "
", predict.cmd, " 
")
    writeLines(script.txt, sh.file)
  }
  cat(sprintf(
    "Writing %4d / %4d samples %d labels %d labeled problems %s\n",
    sample.i, length(sample.dir.vec),
    nrow(labels), length(labels.by.problem),
    problems.dir))
  nothing <- lapply(1:nrow(problems), makeProblem)
  ## Script for peaks on the whole sample.
  peaks.bed <- file.path(sample.dir, "peaks.bed")
  sh.file <- paste0(peaks.bed, ".sh")
  sample.id <- basename(sample.dir)
  script.txt <- paste0(PBS.header, "
#PBS -o ", peaks.bed, ".out
#PBS -e ", peaks.bed, ".err
#PBS -N Predict", sample.id, "
", "Rscript ", normalizePath("predict_sample.R", mustWork=TRUE), " ",
sample.dir, " 
")
  writeLines(script.txt, sh.file)
}

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
cat(
  "Writing ", nrow(problems),
  " jointProblems.bed.sh scripts.",
  "\n", sep="")
for(problem.i in 1:nrow(problems)){
  problem <- problems[problem.i,]
  prob.dir <- file.path(data.dir, "problems", problem$problem.name)
  dir.create(prob.dir, showWarnings=FALSE, recursive=TRUE)
  ## jointProblems.bed.sh
  jointProblems.bed <- file.path(prob.dir, "jointProblems.bed")
  sh.file <- paste0(jointProblems.bed, ".sh")
  pred.cmd <- Rscript(
    'coseg::problem.predict.allSamples("%s")',
    prob.dir)
  joint.prob.cmd <- paste(
    "Rscript",
    normalizePath("create_problems_joint.R", mustWork=TRUE),
    prob.dir)
  joint.targets.cmd <- Rscript(
    'PeakSegJoint::problem.joint.targets("%s")',
    prob.dir)
  script.txt <- paste0(PBS.header, "
#PBS -o ", jointProblems.bed, ".out
#PBS -e ", jointProblems.bed, ".err
#PBS -N P", problem$problem.name, "
", pred.cmd, " 
", joint.prob.cmd, "
", joint.targets.cmd, "
")
  writeLines(script.txt, sh.file)
  ## joint prediction script.
  peaks.bed <- file.path(prob.dir, "peaks.bed")
  sh.file <- paste0(peaks.bed, ".sh")
  pred.cmd <- Rscript(
    'PeakSegJoint::problem.joint.predict.many("%s")',
    prob.dir)
  script.txt <- paste0(PBS.header, "
#PBS -o ", peaks.bed, ".out
#PBS -e ", peaks.bed, ".err
#PBS -N J", problem$problem.name, "
", pred.cmd, " 
")
  writeLines(script.txt, sh.file)
  ## joint prediction jobs script.
  job.dir <- file.path(data.dir, "jobs", problem.i)
  dir.create(job.dir, showWarnings=FALSE, recursive=TRUE)
  job.path <- normalizePath(job.dir, mustWork=TRUE)
  jobPeaks <- file.path(job.dir, "jobPeaks")
  sh.file <- paste0(jobPeaks, ".sh")
  pred.cmd <- Rscript(
    'PeakSegJoint::problem.joint.predict.job("%s")',
    job.path)
  script.txt <- paste0(PBS.header, "
#PBS -o ", jobPeaks, ".out
#PBS -e ", jobPeaks, ".err
#PBS -N Job", problem.i, "
", pred.cmd, " 
")
  writeLines(script.txt, sh.file)
  if(file.exists(chunk.limits.RData) &&
     problem$problem.name %in% chunks.with.problems$problem.name){
    ## jointTargets.sh
    jointTargets <- file.path(prob.dir, "jointTargets")
    sh.file <- paste0(jointTargets, ".sh")
    script.txt <- paste0(PBS.header, "
#PBS -o ", jointTargets, ".out
#PBS -e ", jointTargets, ".err
#PBS -N JT", problem$problem.name, "
", joint.targets.cmd, "
")
    writeLines(script.txt, sh.file)
    ## write a directory for every chunk.
    problem.chunks <- chunks.with.problems[problem$problem.name]
    cat("Writing ", nrow(problem.chunks),
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
      figure.png <- file.path(chunk.dir, "figure-predictions.png")
      script.txt <- paste0(PBS.header, "
#PBS -o ", figure.png, ".out
#PBS -e ", figure.png, ".err
#PBS -N FIG", chunk$chunk.name, "
Rscript ",
normalizePath("plot_chunk.R", mustWork=TRUE),
" ", chunk.dir, "
")
      sh.file <- paste0(figure.png, ".sh")
      writeLines(script.txt, sh.file)
    }
  }
}

## Create joint model script.
joint.model.RData <- file.path(data.dir, "joint.model.RData")
sh.file <- paste0(joint.model.RData, ".sh")
train.cmd <- Rscript(
  'PeakSegJoint::problem.joint.targets.train("%s")',
  data.dir)
script.txt <- paste0(PBS.header, "
#PBS -o ", joint.model.RData, ".out
#PBS -e ", joint.model.RData, ".err
#PBS -N JModel
", train.cmd, " 
")
writeLines(script.txt, sh.file)

## Create plotting script.
peaks.tsv <- file.path(data.dir, "peaks_matrix.tsv")
sh.file <- paste0(peaks.tsv, ".sh")
script.txt <- paste0(
  PBS.header, "
#PBS -o ", peaks.tsv, ".out
#PBS -e ", peaks.tsv, ".err
#PBS -N PeaksMatrix
Rscript ",
  normalizePath("plot_all.R", mustWork=TRUE),
  " ", data.dir, " 
")
writeLines(script.txt, sh.file)

