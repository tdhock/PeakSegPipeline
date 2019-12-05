\name{problem.joint.predict}
\alias{problem.joint.predict}
\title{problem joint predict}
\description{Compute peak predictions for a joint problem.}
\usage{problem.joint.predict(jointProblem.dir, 
    verbose = getOption("PeakSegPipeline.verbose", 
        1))}
\arguments{
  \item{jointProblem.dir}{project/problems/problemID/jointProblems/jointProbID}
  \item{verbose}{
}
}

\value{data.table with one row, describing predicted peaks for both group
and sample models. Columns are chrom, peakStart, peakEnd,
peak.mean.vec, background.mean.vec (for computing normalized peak
height relative to background), sample.loss.diff.vec (for ranking
peaks in a given sample from most to least likely), sample.peaks,
group.peaks, sample.loss.diff, group.loss.diff (for ranking
likelihood of peak regions in joint model).}

\author{Toby Dylan Hocking}





