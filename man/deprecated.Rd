\name{pepStat-deprecated}
\alias{pepStat-deprecated}
\alias{NormalizeArray-deprecated}
\alias{NormalizeArray}
\title{
Deprecated Functions in Package pepStat
}
\description{
The NormalizeArray function was renamed normalizeArray to fall in line
with naming conventions. This function remains deprecated to retain 
functionality in previously written scripts.
}
\usage{
NormalizeArray(peptideSet, robust = TRUE, standard = FALSE, 
  method = "ZpepQuad", centered = TRUE)
}

\arguments{
  \item{peptideSet}{This object contains an ExpressionSet and has the following additional slots}  
  \item{robust}{A \code{logical} value. If TRUE, reweighted least-squares estimates are computed.}
  \item{standard}{A \code{logical}.}
  \item{method}{The normalization method to be used. See Details.}
  \item{centered}{A \code{logical}. If TRUE, recenter the data.}
}

\details{
This function is the same as \link{normalizeArray}, which is what should
be used instead. 
}
\value{
peptideSet with normalized intensities
}
\author{
Gregory Imholte
}

\seealso{
\link{normalizeArray}
}