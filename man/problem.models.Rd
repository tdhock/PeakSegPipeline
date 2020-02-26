\name{problem.models}
\alias{problem.models}
\title{problem models}
\description{Read models from filesystem. Reads both CSV and RDS files, then
saves them all in RDS, then deletes all CSV files.}
\usage{problem.models(problem.dir, 
    verbose = getOption("PeakSegPipeline.verbose", 
        1))}
\arguments{
  \item{problem.dir}{character string path: samples/groupID/sampleID/problems/probID}
  \item{verbose}{print messages?}
}

\value{Data table with all info from model files. One row per model, with
list columns segments.dt and errors.dt with elements that are data
tables representing the segmentation model and label errors.}

\author{Toby Dylan Hocking}





