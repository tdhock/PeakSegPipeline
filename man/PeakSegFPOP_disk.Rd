\name{PeakSegFPOP_disk}
\alias{PeakSegFPOP_disk}
\title{PeakSegFPOP disk}
\description{Run the PeakSeg Functional Pruning Optimal Partitioning algorithm,
using a file on disk (rather than in memory as in
coseg::PeakSegFPOP) to store the O(N) function piece lists, each of size O(log N).
Finds the optimal change-points using the Poisson loss and the
PeakSeg constraint. For N data points, the functional pruning
algorithm is O(N log N) time and disk space, and O(log N) memory.
It computes the exact
solution to the following optimization problem. Let Z be an
N-vector of count data (count.vec, non-negative integers), let W
be an N-vector of positive weights (weight.vec), and let penalty
be a non-negative real number. Find the N-vector M of real numbers
(segment means) and (N-1)-vector C of change-point indicators in
{-1,0,1} which minimize the penalized Poisson Loss,
penalty*sum_{i=1}^{N_1} I(c_i=1) + sum_{i=1}^N
w_i*[m_i-z_i*log(m_i)], subject to constraints: (1) the first
change is up and the next change is down, etc (sum_{i=1}^t c_i in
{0,1} for all t<N-1), and (2) the last change is down
0=sum_{i=1}^{N-1}c_i, and (3) Every zero-valued change-point
variable has an equal segment mean after: c_i=0 implies
m_i=m_{i+1}, (4) every positive-valued change-point variable may
have an up change after: c_i=1 implies m_i<=m_{i+1}, (5) every
negative-valued change-point variable may have a down change
after: c_i=-1 implies m_i>=m_{i+1}. Note that when the equality
constraints are active for non-zero change-point variables, the
recovered model is not feasible for the strict inequality
constraints of the PeakSeg problem, and the optimum of the PeakSeg
problem is undefined.}
\usage{PeakSegFPOP_disk(bedGraph.file, pen.str)}
\arguments{
  \item{bedGraph.file}{character scalar: tab-delimited tabular text file with four
columns: chrom, chromStart, chromEnd, coverage.}
  \item{pen.str}{character scalar that can be converted to a numeric scalar via
as.numeric: non-negative penalty. More penalty means fewer
peaks. 0 and Inf are OK.}
}



\author{Toby Dylan Hocking}




\examples{

library(PeakSegPipeline)
r <- function(chrom, chromStart, chromEnd, coverage){
  data.frame(chrom, chromStart, chromEnd, coverage)
}
four <- rbind(
  r("chr1", 0, 10,  2),
  r("chr1", 10, 20, 10),
  r("chr1", 20, 30, 14),
  r("chr1", 30, 40, 13))
write.table(
  four, tmp <- tempfile(),
  sep="\t", row.names=FALSE, col.names=FALSE)
names.list <- PeakSegFPOP_disk(tmp, "10.5")
unlink(names.list$db)
seg.df <- read.table(names.list$segments)
names(seg.df) <- c("chrom", "chromStart", "chromEnd", "status", "mean")
seg.df
loss.df <- read.table(names.list$loss)
names(loss.df) <- c(
  "penalty", "segments", "peaks", "bases",
  "mean.pen.cost", "total.cost", "status",
  "mean.intervals", "max.intervals")
loss.df

}
