\name{problem.PeakSegFPOP}
\alias{problem.PeakSegFPOP}
\title{problem PeakSegFPOP}
\description{Run PeakSegFPOP on one genomic segmentation problem directory.}
\usage{problem.PeakSegFPOP(problem.dir, penalty.str)}
\arguments{
  \item{problem.dir}{Path to a directory like sampleID/problems/problemID which
contains a coverage.bedGraph file with the aligned read counts for
one genomic segmentation problem.}
  \item{penalty.str}{Penalty parameter to pass to the PeakSegFPOP command line program.}
}

\value{List of data.tables: segments has one row for every segment in the
optimal model, loss has one row and contains the Poisson loss and
feasibility, and timing is one row with the time and disk usage.}

\author{Toby Dylan Hocking}





