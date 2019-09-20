\name{downloadBigWigs}
\alias{downloadBigWigs}
\title{downloadBigWigs}
\description{Download bigWig files from a trackDb file.}
\usage{downloadBigWigs(trackDb.txt, out.dir, verbose = 0)}
\arguments{
  \item{trackDb.txt}{trackDb text file.}
  \item{out.dir}{Output directory.}
  \item{verbose}{Print messages?}
}

\value{Nothing.}

\author{Toby Dylan Hocking}




\examples{

## Create a trackDb.txt file with links to a few small bigWig files on GitHub.
bigWig.part.vec <- c(
  "Input/MS010302",
  "bcell/MS010302",
  "Input/MS002201",
  "kidney/MS002201")
repos.url <- "https://raw.githubusercontent.com/tdhock/input-test-data/master/"
track.lines <- sprintf(
  "track \%s\\nbigDataUrl \%s\%s.bigwig\\n",
  sub("/", "_", bigWig.part.vec),
  repos.url, bigWig.part.vec)
track.dir <- tempfile()
dir.create(track.dir)
trackDb.txt <- file.path(track.dir, "trackDb.txt")
cat(track.lines, sep="\\n")
writeLines(track.lines, trackDb.txt)

## Download the bigWig files mentioned in trackDb.txt
if(interactive()){
  downloadBigWigs(trackDb.txt, track.dir, verbose=1)
}

}
