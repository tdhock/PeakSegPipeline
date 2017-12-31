\name{create_track_hub}
\alias{create_track_hub}
\title{create track hub}
\description{Create track hub for a project.}
\usage{create_track_hub(data.dir.path, data.dir.url, genome, 
    email, goldenPath.url = ucsc.goldenPath.url)}
\arguments{
  \item{data.dir.path}{data/project directory path on local filesystem.}
  \item{data.dir.url}{URL for data/project directory, for constructing URLs to
bigWig/bigBed files in the trackDb.txt file, e.g. if
data.dir.url="http://some.domain/~user/demo_project" then your
trackDb.txt file will contain URLs like
http://some.domain/~user/demo_project/samples/groupID/sampleID/coverage.bigWig}
  \item{genome}{genome string as defined at UCSC, e.g. "hg19"}
  \item{email}{email address for maintainer of track hub.}
  \item{goldenPath.url}{link to download UCSC genome chromInfo files, necessary for
creating bigWigs.}
}



\author{Toby Dylan Hocking}





