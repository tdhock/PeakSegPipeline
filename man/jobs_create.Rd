\name{jobs_create}
\alias{jobs_create}
\title{jobs create}
\description{Setup a data directory for analysis with PeakSegPipeline.}
\usage{jobs_create(data.dir.arg, 
    target.minutes = NULL, 
    verbose = getOption("PeakSegPipeline.verbose", 
        1))}
\arguments{
  \item{data.dir.arg}{path to project directory.}
  \item{target.minutes}{Maximum number of minutes to spend computing the target interval
during Step1 of the pipeline (learning what penalties result in
low label error). If this is a number then a \code{target.minutes} file
will be written to each problem directory -- this number in this
file is used as the maximum number of minutes by the
\code{\link{problem.target}} function.}
  \item{verbose}{TRUE for output, FALSE otherwise.}
}

\value{data.table with one row for each job and three columns: fun, arg,
step. fun is the function to call with argument arg, in order
specified by step (smaller steps first).}

\author{Toby Dylan Hocking}





