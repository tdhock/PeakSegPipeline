\name{create_track_hub}
\alias{create_track_hub}
\title{create track hub}
\description{Create track hub for a project.}
\usage{create_track_hub(data.dir = "test/input", url.prefix = "http://hubs.hpc.mcgill.ca/~thocking/PeakSegFPOP-", 
    genome = "hg19", email = "email@domain.com", goldenPath.url = "http://hgdownload.soe.ucsc.edu/goldenPath/")}
\arguments{
  \item{data.dir}{data/project directory.}
  \item{url.prefix}{Prefix to use for links to bigWig/bigBed files, data.dir will be
appended after this. e.g. if
url.prefix="http://some.domain/~user/foo-" and
data.dir="test/input" then URLS will be
http://some.domain/~user/foo-test/input/samples/groupID/sampleID/coverage.bigWig,
etc.}
  \item{genome}{genome string as defined at UCSC, e.g. "hg19"}
  \item{email}{email address for maintainer of track hub.}
  \item{goldenPath.url}{link to download UCSC genome chromInfo files, necessary for
creating bigWigs.}
}



\author{Toby Dylan Hocking}





