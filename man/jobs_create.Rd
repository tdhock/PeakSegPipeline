\name{jobs_create}
\alias{jobs_create}
\title{jobs create}
\description{Setup a data directory for analysis with PeakSegPipeline.}
\usage{jobs_create(data.dir.arg, 
    verbose = getOption("PeakSegPipeline.verbose", 
        1))}
\arguments{
  \item{data.dir.arg}{path to project directory.}
  \item{verbose}{TRUE for output, FALSE otherwise.}
}

\value{data.table with one row for each job and three columns: fun, arg,
step. fun is the function to call with argument arg, in order
specified by step (smaller steps first).}

\author{Toby Dylan Hocking}





