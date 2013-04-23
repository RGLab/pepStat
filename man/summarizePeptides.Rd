\name{summarizePeptides}
\alias{summarizePeptides}

\title{Add informations to a peptideSet}

\description{
This function merges the replicates and adds information from a peptide collection to a peptideSet. This include coordinates, alignment information and zScores.
}

\usage{
summarizePeptides(peptideSet,summary="median",position=NULL)
}

\arguments{
\item{peptideSet}{A \code{peptideSet} object.}
\item{summary}{A character string. The method to be used for merging replicates.}
\item{position}{A \code{RangedData} object. A collection of peptides that can be found in the PEP.db package.}
%\item{\dots}{Additional items.}
}

\value{
An object of class \code{peptideSet} with added columns and updated ranges.
}

\note{
For the position argument, some peptide collections can be found in the \code{PEP.db} package.
}

\author{Raphael Gottardo}

\seealso{
\code{makePeptideSet}
}

\keyword{functions}
