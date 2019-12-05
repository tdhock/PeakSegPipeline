\name{problem.joint.predict.many}
\alias{problem.joint.predict.many}
\title{problem joint predict many}
\description{Compute all joint peak predictions for one separate problem, in
parallel over joint problems using \code{\link[future.apply]{future_lapply}}.}
\usage{problem.joint.predict.many(prob.dir, 
    verbose = getOption("PeakSegPipeline.verbose", 
        1))}
\arguments{
  \item{prob.dir}{project/problems/problemID}
  \item{verbose}{
}
}

\value{data.table of predicted peaks.}

\author{Toby Dylan Hocking}





