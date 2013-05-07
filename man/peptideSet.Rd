\name{peptideSet-class}
\docType{class}
\alias{peptideSet}
\alias{peptideSet-class}
\alias{clade,peptideSet-method}
\alias{[,peptideSet,ANY,ANY,logical-method}
\alias{subset,peptideSet-method}
\alias{start,peptideSet-method}
\alias{end,peptideSet-method}
\alias{width,peptideSet-method}
\alias{values,peptideSet-method}
\alias{values<-,peptideSet-method}
\alias{show,peptideSet-method}
\alias{summary,peptideSet-method}
\alias{peptide}
\alias{peptide,peptideSet-method}
\alias{position}
\alias{position,peptideSet-method}
\alias{split,peptideSet-method}
\alias{featureID}
\alias{featureID,peptideSet-method}
\alias{write.pSet}
\alias{write.pSet,peptideSet-method}


\title{The peptideSet class}

\description{
This object gathers all information from gpr files, annotation data, and sequence data.
}

\section{Accessors}{
The pepStat package offers accessors to most of the attributes of a peptideSet object:
\describe{
\item{start(x)}{Returns the start of each peptide.}
\item{end(x)}{Returns the end of each peptide.}
\item{width(x)}{Returns the length of each peptide.}
\item{position(x)}{Returns the coordinates of the central amino-acid of each peptide, equal to \code{round((start(x) + end(x))/2)}}
\item{ranges(x)}{Returns a \code{RangedData} object.}
\item{values(x)}{Returns a \code{CompressedSplitDataFrameList} object.}
}
}

\section{Constructor}{
The \link{makePeptideSet} function can be used to create a peptideSet object from a list of files or a folder.
}

\section{Methods}{
\describe{
  \item{"["}{: \code{signature(x = "peptideSet")}: subset method for \code{peptideSet} object.}
  \item{"subset"}{: \code{signature(x = "peptideSet")}: subset method for \code{peptideSet} object.}
  \item{"clade"}{: \code{signature(object = "peptideSet")}: Returns a \code{matrix} of the clades each peptide belongs to.}
  \item{"show"}{: \code{signature(x = "peptideSet")}: Display a quick summary of the Samples and the preprocessing information.}
  \item{"peptide"}{: \code{signature(x = "peptideSet")}: Returns the sequence of the peptides.}
  \item{"split"}{: \code{signature(x = "peptideSet", f = "factor")}: Split the peptideSet into a \code{list} of \code{peptideSet} based on the given factor.}
  \item{"featureID"}{: \code{signature(x = "peptideSet")}: The ID of the peptides.}
  \item{"write.pSet"}{: \code{signature(x = "peptideSet")}: Write the peptides to a file. This is a wrapper to \code{write.table} and any additional argument \code{"..."} will be passed to the base function.}
}
}




\author{Raphael Gottardo}


\seealso{
\code{\link{makePeptideSet}}
}

\keyword{models}
