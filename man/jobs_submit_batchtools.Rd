\name{jobs_submit_batchtools}
\alias{jobs_submit_batchtools}
\title{jobs submit batchtools}
\description{Submit PeakSegPipeline jobs via batchtools.}
\usage{jobs_submit_batchtools(jobs, resources = list(walltime = 24 * 
    60, memory = 2000, ncpus = 2, ntasks = 1, chunks.as.arrayjobs = TRUE))}
\arguments{
  \item{jobs}{data.table from jobs_create.}
  \item{resources}{
}
}



\author{Toby Dylan Hocking}




\examples{
if(FALSE){
  jobs <- jobs_create("~/genomic-ml/PeakSegFPOP/labels/ATAC_JV_adipose/")
  jobs_submit_batchtools(jobs)
}
}