\name{PPN.cores}
\alias{PPN.cores}
\title{PPN cores}
\description{Set mc.cores option from the PBS_NUM_PPN environment variable.}
\usage{PPN.cores(variable = "PBS_NUM_PPN")}
\arguments{
  \item{variable}{
}
}

\value{integer: ppn value obtained from environment variable -- it will
only be used to set the mc.cores option if it is finite. (if the
env var is undefined then this will be NA so this function will
not set the mc.cores option)}

\author{Toby Dylan Hocking}





