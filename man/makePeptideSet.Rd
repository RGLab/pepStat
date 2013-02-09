\name{makePeptideSet}
\alias{makePeptideSet}

\title{peptideSet constructor}

\description{
This function reads raw microarray data files and create a peptideSet object to gather the information.
}

\usage{
makePeptideSet(files=NULL, path=NULL, mapping.file=NULL, use.flags=FALSE, rm.control.list=c("empty","none","JPT-","Ig","Cy","landmark"), norm.empty=TRUE, empty.control.list=c("empty","blank control"), bgCorrect.method="normexp", log=TRUE, check.row.order=FALSE, verbose=FALSE)
}

\arguments{
\item{files}{A character vector containing the names of the files to read. If omitted, all the files with the extension .gpr in specified folder will be read.}
\item{path}{A character string giving the directory containing files.}
\item{mapping.file}{A character string. The name of a mapping file for the data. This will add the relevant phenoData to the peptideSet.}
\item{use.flags}{A \code{logical}.}
\item{rm.control.list}{A vector of character strings. The name of the controls to be excluded from the peptideSet.}
\item{norm.empty}{A \code{logical}.}
\item{empty.control.list}{A vector of character strings. The name of the empty controls. These are used or background correction.}
\item{bgCorrect.method}{A character string. The name of the method used for background correction.}
\item{log}{A \code{logical}.}
\item{check.row.order}{A \code{logical}.}
\item{verbose}{A \code{logical}.}
}

\author{Raphael Gottardo}

\seealso{
\code{peptideSet}
}

\keyword{functions}
