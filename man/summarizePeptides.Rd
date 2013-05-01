\name{summarizePeptides}
\alias{summarizePeptides}

\title{Add informations to a peptideSet and summarizes peptides}

\description{
This function merges the replicates and adds information from a peptide collection 
to a peptideSet. This include coordinates, alignment information and zScores.
}

\usage{
summarizePeptides(peptideSet,summary="median",position=NULL)
}

\arguments{
\item{peptideSet}{A \code{peptideSet} object, preferably created by \link{makePeptideSet}.}
\item{summary}{A character string. The method to be used for merging replicates.}
\item{position}{A \code{RangedData} object. See Details}
}

\details{
The \code{summary} argument must be either "mean" or "median", indicating how within-array
replicates of peptides will be summarized.

If non-\code{NULL}, the \code{position} argument must be a \link{RangedData} object 
that has rownames corresponding to the ID of array probes. In particular,
array probe IDs from \code{peptide(peptideSet)} are matched against these rownames. 
It is permissible for the \code{position} object to cover multiple spaces.
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
