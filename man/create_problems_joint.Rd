\name{create_problems_joint}
\alias{create_problems_joint}
\title{create problems joint}
\description{Create joint problems for one separate problem, after separate
peak prediction. Parallelized over joint problems via
\code{\link{psp_lapply}}.}
\usage{create_problems_joint(prob.dir, 
    peaks = NULL, verbose = getOption("PeakSegPipeline.verbose", 
        1))}
\arguments{
  \item{prob.dir}{proj.dir/problems/problemID}
  \item{peaks}{data.table of \code{peaks} predicted in all samples for this problem (if
it has already been computed by \code{\link{problem.predict.allSamples}}), or
NULL which means to read predicted \code{peaks} from
proj.dir/samples/*/*/problemID/peaks.bed files.}
  \item{verbose}{
}
}



\author{Toby Dylan Hocking}





