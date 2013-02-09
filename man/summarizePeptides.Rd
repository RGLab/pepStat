\name{summarizePeptides}
\alias{summarizePeptides}

\title{Add informations to a peptideSet}

\description{
This function merges the replicates and adds information from a peptide collection to a peptideSet. This include coordinates, alignment information and zScores.
}

\usage{
summarizePeptides<-function(peptideSet,summary="mean",position=NULL,...)
}

\arguments{
\item{peptideSet}{A \code{peptideSet} object.}
\item{summary}{A character string. The method to be used for merging replicates.}
\item{position}{A \code{RangedData} object. A collection of peptides that can be found in the PEP.db package.}
}

\value{
An object of class \code{peptideSet} with added columns and updated ranges.
}

\author{Raphael Gottardo}

\seealso{
\code{makePeptideSet}
}

\keyword{functions}
