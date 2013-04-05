\name{peptideSet-class}
\docType{class}
\alias{peptideSet}
\alias{peptideSet-class}
\alias{clade,peptideSet-method}
\alias{[,peptideSet-method}
\alias{subset,peptideSet-method}
\alias{start,peptideSet-method}
\alias{end,peptideSet-method}
\alias{width,peptideSet-method}
\alias{values,peptideSet-method}
\alias{values<-,peptideSet-method}
\alias{show,peptideSet-method}
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
This object gather all information from gpr files.
}

\section{Accessors}{
The pepStat package offers accessors to most of the attributes of a peptideSet object:
\describe{
\item{start(x)}{Returns the start of each peptide.}
\item{end(x)}{Returns the end of each peptide.}
\item{width(x)}{Returns the length of each peptide.}
\item{position(x)}{Returns the coordinates of the central amino-acid of each peptide.}
\item{ranges(x)}{Returns a \code{RangedData} object.}
\item{values(x)}{Returns a \code{CompressedSplitDataFrameList} object.}
}
}

\section{Constructor}{
The makePeptideSet function can be used to create a peptideSet object from a list of file or a folder.
}

\section{Methods}{
\describe{
  \item{"["}{: \code{signature(x = "peptideSet")}: subset method for \code{peptideSet} object.}
  \item{"subset"}{: \code{signature(x = "peptideSet")}: subset method for \code{peptideSet} object.}
  \item{"clade"}{: \code{signature(object = "peptideSet")}: Returns a \code{matrix} of the clades each peptide belongs to.}
  %\item{"start"}{: \code{signature(x = "peptideSet")}: Returns the starting position of the peptides.}
  %\item{"end"}{: \code{signature(x = "peptideSet")}: Returns the last position of the peptides.}
  %\item{"width"}{: \code{signature(x = "peptideSet")}: Returns the length of the peptides.}
  %\item{"values"}{:}
  %\item{"values<-"}{:}
  \item{"show"}{: \code{signature(x = "peptideSet")}: Display method.}
  \item{"peptide"}{: \code{signature(x = "peptideSet")}: Returns the sequence of the peptides.}
  %\item{"position"}{: \code{signature(x = "peptideSet")}: The position of the peptides: (end - start)/2 }
  \item{"split"}{: \code{signature(x = "peptideSet", f = "factor")}: Split the peptideSet in a \code{list} of \code{peptideSet} based on the given factor.}
  \item{"featureID"}{: \code{signature(x = "peptideSet")}: The ID of the peptides.}
  \item{"write.pSet"}{: \code{signature(x = "peptideSet")}: Write the peptides to a file. This is a wrapper to \code{write.table} and any additional argument \code{"..."} will be passed to the base function.}
}
}




\author{Raphael Gottardo}


\seealso{
\code{\link{makePeptideSet}}
}

\keyword{models}
