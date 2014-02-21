\name{summarizePeptides}
\alias{summarizePeptides}
\title{Add information to a peptideSet and summarize peptides}
\usage{
summarizePeptides(peptideSet, summary = "median", position = NULL)
}
\arguments{
  \item{peptideSet}{A \code{peptideSet}, as created by
  \code{makePeptideSet}}

  \item{summary}{A \code{character} string. The method used
  for merging replicates. Available are: "mean" and
  "median".}

  \item{position}{A \code{RangedData}. A peptide collection
  such as the ones available in \code{PEP.db}. See details
  below and vignettes for more information.}
}
\value{
An object of class \code{peptideSet} with added columns and
updated ranges.
}
\description{
This function merges the replicates and adds information
from a peptide collection to a peptideSet. This collection
can include coordinates, alignment information, Z-scales,
and other peptide information.
}
\note{
For the position argument, some peptide collections can be
found in the \code{PEP.db} package.
}
\author{
Raphael Gottardo, Greory Imholte
}
\seealso{
\code{\link{makePeptideSet}}, \code{\link{create_db}}
}

