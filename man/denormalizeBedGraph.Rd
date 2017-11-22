\name{denormalizeBedGraph}
\alias{denormalizeBedGraph}
\title{denormalizeBedGraph}
\description{Attempt to convert a bedGraph file with non-negative non-integer
data to integer count data, by dividing each value by the smallest
non-zero value.}
\usage{denormalizeBedGraph(input.bedGraph, output.bedGraph)}
\arguments{
  \item{input.bedGraph}{Input bedGraph file with non-integer data.}
  \item{output.bedGraph}{Output bedGraph file with integer data.}
}

\value{Nothing.}

\author{Toby Dylan Hocking}





