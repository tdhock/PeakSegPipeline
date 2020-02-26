\name{problem.labels}
\alias{problem.labels}
\title{problem labels}
\description{read problemID/labels.bed if it exists, otherwise read
sampleID/labels.bed}
\usage{problem.labels(problem.dir, 
    stop.without.labels = TRUE)}
\arguments{
  \item{problem.dir}{project/samples/groupID/sampleID/problems/problemID}
  \item{stop.without.labels}{
}
}

\value{data.table with one row for each label and columns chrom,
chromStart, chromEnd, annotation.}

\author{Toby Dylan Hocking}





