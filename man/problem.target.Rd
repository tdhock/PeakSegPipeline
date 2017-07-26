\name{problem.target}
\alias{problem.target}
\title{problem target}
\description{Compute target interval for a segmentation problem.}
\usage{problem.target(problem.dir)}
\arguments{
  \item{problem.dir}{problemID directory in which coverage.bedGraph has already been
computed. If there is a labels.bed file then the number of
incorrect labels will be computed in order to find the target
interval of minimal error penalty values.}
}

\value{List of info related to target interval computation: target is the
interval of log(penalty) values that achieve minimum incorrect
labels (numeric vector of length 2), models is a data.table with
one row per model for which the label error was computed, selected
is a data.frame that describes which penalty values will select
which models.}

\author{Toby Dylan Hocking}





