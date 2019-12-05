\name{problem.joint.targets}
\alias{problem.joint.targets}
\title{problem joint targets}
\description{Compute joint targets for all joint problems in a separate
problem, in parallel over joint problems using
\code{\link[future.apply]{future_lapply}}.}
\usage{problem.joint.targets(problem.dir, 
    verbose = getOption("PeakSegPipeline.verbose", 
        1))}
\arguments{
  \item{problem.dir}{project/problems/problemID}
  \item{verbose}{
}
}

\value{Named list of two matrices: targets and features.}

\author{Toby Dylan Hocking}





