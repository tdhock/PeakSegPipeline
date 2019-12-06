\name{gap2problems}
\alias{gap2problems}
\title{gap2problems}
\description{Convert a gap file from UCSC to a \code{problems.bed} file, which is used
in the PeakSegPipeline to determine where PeakSegFPOP should be
run to do peak calling for each sample independently.}
\usage{gap2problems(gap.bed, 
    chromInfo.txt, problems.bed)}
\arguments{
  \item{gap.bed}{\code{gap.bed} file, parts of the genome with gaps (NNNN).}
  \item{chromInfo.txt}{\code{chromInfo.txt} file, chromosome sizes.}
  \item{problems.bed}{will be created}
}



\author{Toby Dylan Hocking}





