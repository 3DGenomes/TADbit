\name{tadbit}
\alias{tadbit}

\title{Identify TADs by breakpoint detection}
\description{
Identify topologically associated domains in hi-C data by optimal 
breakpoint detection.
}
\usage{
tadbit(x, max_size="auto", n_CPU="auto", verbose=FALSE)
}

\arguments{
  \item{x}{A square \code{matrix} of interaction counts in hi-C data or
  a \code{list} of such matrices for replicated experiments. The counts
  must be evenly sampled and not normalized.}
  \item{max_size}{The size of the largest possible TAD expressed either
  as a fraction of the length of the chromosome, or in absolute number
  of count bins. The value "auto" is 10\% of the chromosome length.}
  \item{n_CPU}{The number of CPUs to allocate to tadbit. The value "auto"
  is the total number of CPUs minus 1.}
  \item{verbose}{Whether to indicate progression of the algorithm.}
}

\details{
  The tadbit algorithm works on raw chromosome interaction count data.
  Not only is normalization not necessary, it is also not recommended
  since the data is assumed to be discrete counts.

  Tadbit is a breakpoint detection algorithm that returns the optimal
  segmentation of the chromosome under BIC-penalized likelihood. The
  model assumes that counts have a Poisson distribution and that the
  expected value of the counts decreases like a power-law with the
  linear distance on the chromosome. This expected value of the counts
  at position (i,j) is corrected by the counts at diagonal positions
  (i,i) and (j,j). This normalizes for different restriction enzynme
  site densities and 'mappability' of the reads in case a bin contains
  repeated regions.
}

\value{
  A \code{vector} indicating the breaks in the original data.
}

\author{Guillaume Filion}
\keyword{breakpoint detection}
\keyword{hi-C}

