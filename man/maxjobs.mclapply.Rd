\name{maxjobs.mclapply}
\alias{maxjobs.mclapply}
\title{maxjobs mclapply}
\description{Run mclapply inside of a for loop, ensuring that it never receives
a first argument with a length more than maxjobs. This avoids some
memory problems (swapping, or getting jobs killed on the cluster)
when using mclapply(1:N, FUN) where N is large.}
\usage{maxjobs.mclapply(X, FUN, maxjobs = getOption("mc.cores", 
    1L))}
\arguments{
  \item{X}{
}
  \item{FUN}{
}
  \item{maxjobs}{
}
}



\author{Toby Dylan Hocking}





