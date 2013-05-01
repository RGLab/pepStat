\name{NormalizeArray}
\alias{NormalizeArray}

\title{Normalize tiling array data using sequence information}

\description{
  This function is used to normalize peptide microarray data using sequence information. Users can chose between four different normalization methods. Please refer to the arguments section below. 
}

\usage{
NormalizeArray(peptideSet, robust = TRUE, standard = FALSE, method = "ZpepQuad",
  centered = TRUE, verbose = FALSE)
}

\arguments{
  \item{peptideSet}{This object contains an ExpressionSet and has the following additional slots}	
  \item{robust}{A \code{logical} value. If TRUE, reweighted least-squares estimates are computed.}
  \item{standard}{A \code{logical}.}
  \item{method}{The normalization method to be used. See Details.}
  \item{centered}{A \code{logical}. If TRUE, recenter the data.}
  \item{verbose}{A \code{logical} value. If verbose is TRUE, progress information is be displayed.}
}

\details{

}

\value{
  A \link{peptideSet} object with updated normalized intensity values.
}

\author{
  Raphael Gottardo, Gregory Imholte
}
