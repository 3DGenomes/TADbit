\name{batch_tadbit}
\alias{batch_tadbit}

\title{Batch mode tadbit}
\description{
Use tadbit on directories of data files
}
\usage{
batch_tadbit(directory=getwd(), sep='_', read_options=list(), ...)
}

\arguments{
  \item{directory}{The directory containing the data files.}
  \item{sep}{A character specifying how to identify unit/chormosome names
  (see below).}
  \item{read_options}{A \code{list} of options passed to \code{read.delim}
  upon reading data files.}
  \item{...}{Further arguments passed to \code{tadbit}.}
}

\details{
  All files in the specified directory will be considered data file. The
  presence of non data files will cause the function to either crash or
  produce aberrant results.

  Each file has to contain the data for a single unit/chromosome. The
  files can be separated in sub-directories corresponding to single
  experiments or any other organization. Data files that should be
  considered replicates have to start with the same characters, until
  the character \code{sep}. For instance, all replicates of the unit
  "chr1" should start with "chr1_", using the default value of \code{sep}.

  The data files are read through \code{read.delim}. You can pass options
  to \code{read.delim} through the list \code{read_options}. For instance
  if the files have no header, use \code{read_options=list(header=FALSE)} and if
  they also have row names, \code{read_options=list(header=FALSE, row.names=1)}.

  Other arguments such as \code{max_size}, \code{n_CPU} and \code{verbose} are
  passed to \code{tadbit}.
}

\value{
  A \code{list} where each element has the name of the unit/chromosome,
  and is the output of \code{tadbit} run on the corresponding files
  assumed to be replicates.
}

\seealso{
  \code{tadbit}
}

\author{Guillaume Filion}
\keyword{breakpoint detection}
\keyword{hi-C}

