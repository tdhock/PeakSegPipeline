\name{downloadBigWigs}
\alias{downloadBigWigs}
\title{downloadBigWigs}
\description{Download bigWig files from a trackDb file.}
\usage{downloadBigWigs(trackDb.txt, out.dir)}
\arguments{
  \item{trackDb.txt}{trackDb text file.}
  \item{out.dir}{Output directory.}
}

\value{Nothing.}

\author{Toby Dylan Hocking}




\examples{

trackDb.txt <- paste0(
  "https://rcdata.nau.edu",
  "/genomic-ml/PeakSegFPOP/labels/H3K36me3_TDH_immune/trackDb.txt")
downloadBigWigs(trackDb.txt, tempfile())

}
