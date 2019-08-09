\name{plot_all}
\alias{plot_all}
\title{plot all}
\description{Gather and plot results of peak calling, generate summary web page
set.dir.arg/index.html. Labeled chunk plots are created in
parallel via \code{\link[future.apply]{future_lapply}}. If set.dir.arg/hub.sh
exists it is called at the end of this function in order to
generate a track hub based on the peak model files -- it should
contain something like \code{\link{Rscript}} -e
'PeakSegPipeline::create_track_hub(...)'}
\usage{plot_all(set.dir.arg, zoom.out.times = 10)}
\arguments{
  \item{set.dir.arg}{Path/to/data/dir.}
  \item{zoom.out.times}{
}
}

\value{Nothing.}

\author{Toby Dylan Hocking}





