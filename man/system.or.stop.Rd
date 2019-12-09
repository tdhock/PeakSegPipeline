\name{system.or.stop}
\alias{system.or.stop}
\title{system or stop}
\description{Run a command line and stop with
an error for non-zero status code.}
\usage{system.or.stop(cmd, verbose = getOption("PeakSegPipeline.verbose", 
    1))}
\arguments{
  \item{cmd}{Command line arguments, which will be escaped by \code{\link[base]{shQuote}},
pasted together, and then passed to system.}
  \item{verbose}{print output?}
}

\value{Nothing.}

\author{Toby Dylan Hocking}





