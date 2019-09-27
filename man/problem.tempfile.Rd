\name{problem.tempfile}
\alias{problem.tempfile}
\title{problem tempfile}
\description{Create a (problem,penalty)-specific temporary file name to pass to
PeakSegFPOP_dir as the cost function database. The
PeakSegPipeline.TMPDIR option is used for the base temporary
directory if that is set, and otherwise tempdir() is used.}
\usage{problem.tempfile(problem.dir, pen.str)}
\arguments{
  \item{problem.dir}{full path to problem directory.}
  \item{pen.str}{penalty string.}
}

\value{character: temporary file name.}

\author{Toby Dylan Hocking}





