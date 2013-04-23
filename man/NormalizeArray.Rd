\name{NormalizeArray}
\alias{NormalizeArray}

\title{Normalize tiling array data using sequence information}

\description{
  This function is used to normalize tiling array data using sequence information. Users can chose between two different normalization methods. Please refer to the arguments section below. 
}

\usage{
NormalizeArray(peptideSet, robust=TRUE, standard=FALSE, method="ZpepQuad",centered=TRUE, verbose=FALSE)
}
%ScSetNorm<-NormalizeProbes(tilingSet, method="MAT",robust=FALSE, all=FALSE, standard=TRUE, verbose=FALSE)

\arguments{
  %\item{tilingSet}{This object contains an ExpressionSet and has the following additional slots}	
  \item{peptideSet}{This object contains an ExpressionSet and has the following additional slots}	
  \item{robust}{A \code{logical} value. If TRUE, reweighted least-squares estimates are computed.}
  \item{standard}{A \code{logical}.}
  \item{method}{The normalization method to be used. User can choose from "MAT", or "PairBinned". As an upgrade to MAT, the Pair option also takes into account of the interaction between adjacent pairs along the probe as covariates for linear regression. }
  \item{centered}{A \code{logical}. If TRUE, recenter the data.}
  %\item{all}{A logical value. If not using all probes to compute (for faster computation and memory efficiency) the regression parameters, then use the minimum of 300,000 or number of probes, whichever is less.}
  \item{verbose}{A \code{logical} value. If verbose is TRUE, progress information would be displayed.}
}

\details{
For the original rMAT normalization: method is set to be rMAT in string, robust is set to be false, copyNumber is set to be your copy number's vector, rMATScaling is set to be true, and logTransform is set to be true for untransformed data. The output can be saved as BAR file if the BAR argument specifies a filename, or as a parsed BAR file if argument output specifies a filename. 


For more details on normalization, please refer to the following paper:
Johnson et al. Model-based analysis of tiling-arrays for ChIP-chip. Proc Natl Acad Sci USA (2006) vol. 103 (33) pp. 12457-62

}

\value{
  The matrix of normalized expression values.
}

\author{
  Charles Cheung, \email{cykc@interchange.ubc.ca} and
  Raphael Gottardo, \email{raphael.gottardo@ircm.qc.ca}
  Arnaud Droit, \email{arnaud.droit@ircm.qc.ca}
}


\seealso{
  \code{PairInMatrix()} for generating neighbouring pair-codes from sequences and \code{affyTile} for information about the package.
}

\examples{

####################################################
#The data are in inst/doc folder in rMAT package.
####################################################



}

\keyword{file}
\keyword{IO}
