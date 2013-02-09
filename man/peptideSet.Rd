\name{peptideSet-class}
\docType{class}
\alias{peptideSet}
\alias{peptideSet-class}


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

\author{Raphael Gottardo}


\seealso{
\code{\link{makePeptideSet}}
}

\keyword{models}
