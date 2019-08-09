\name{problem.joint.predict.job}
\alias{problem.joint.predict.job}
\title{problem joint predict job}
\description{Compute all joint peak predictions for the joint problems listed
in jobProblems.bed, in parallel over problems using
\code{\link[future.apply]{future_lapply}}.}
\usage{problem.joint.predict.job(job.dir)}
\arguments{
  \item{job.dir}{project/jobs/jobID}
}

\value{data.table of predicted peaks, one row for each job, same columns
as from \code{\link{problem.joint.predict}}.}

\author{Toby Dylan Hocking}





