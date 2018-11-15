\name{jobs_submit_mclapply}
\alias{jobs_submit_mclapply}
\title{jobs submit mclapply}
\description{Run PeakSegPipeline jobs in this R session,
parallelizing the jobs in each step via mclapply.or.stop.}
\usage{jobs_submit_mclapply(jobs)}
\arguments{
  \item{jobs}{data.table from jobs_create.}
}



\author{Toby Dylan Hocking}




\examples{
if(FALSE){
  jobs <- jobs_create("~/genomic-ml/PeakSegFPOP/labels/ATAC_JV_adipose/")
  jobs_submit_mclapply(jobs)
}
}
