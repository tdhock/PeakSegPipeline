\name{create_problems_joint}
\alias{create_problems_joint}
\title{create problems joint}
\description{Create joint problems for one separate problem, after separate
peak prediction.}
\usage{create_problems_joint(prob.dir, peaks = NULL)}
\arguments{
  \item{prob.dir}{proj.dir/problems/problemID}
  \item{peaks}{data.table of peaks predicted in all samples for this problem (if
it has already been computed by problem.predict.allSamples), or
NULL which means to read predicted peaks from
proj.dir/samples/*/*/problemID/peaks.bed files.}
}



\author{Toby Dylan Hocking}





