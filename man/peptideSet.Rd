\name{peptideSet}
\alias{peptideSet}
\alias{peptideSet-class}
\title{peptideSet class}
\description{
  This class gathers all information from gpr files,
  annotation data and sequence data
}
\section{Slots}{
  \describe{ \item{featureRange:}{ A \code{RangedData
  object}. The ranges and sequences of the peptides and
  their associated annotation.} \item{phenoData:}{ An
  \code{AnnotatedDataFrame}. Annotation for the samples.} }
}
\author{
  Greg Imholte
}
\seealso{
  \code{\link{ExpressionSet}},
  \code{\link{peptideSet-methods}}
}

