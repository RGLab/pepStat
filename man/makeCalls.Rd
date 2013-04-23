\name{makeCalls}
\alias{makeCalls}

\title{Making antibody binding positivity calls}

\description{
After normalization and data smoothing, this last step makes the calls for each peptide of the peptideSet.
}

\usage{
makeCalls(peptideSet, cutoff=.1, method="local", freq=TRUE, group=NULL, verbose=FALSE)
}

\arguments{
  \item{peptideSet}{A \code{peptideSet} object. The peptides, after normalization and data smoothing.}
  \item{cutoff}{A \code{numeric}.}
  \item{method}{A \code{character}. The method used to make the positivity calls.}
	\describe{
		\item{local}{}
		\item{global}{}
		\item{absolute}{A simple threshold. All values over \code{cutoff} are identified as binding events.}
		}
  \item{freq}{A \code{logical}. }
  \item{group}{A \code{character}.}
  \item{verbose}{A \code{logical}. If set to TRUE, progress information will be displayed.}
}

\value{
The function returns a \code{numeric} vector with one value for each peptide. A value superior to zero indicates a binding event.
}


\seealso{
\code{NormalizeArray}
\code{slidingMean}
}

\keyword{functions}
