\name{NormalizeArray}
\alias{NormalizeArray}

\title{Normalize tiling array data using sequence information}

\description{
  This function is used to normalize peptide microarray data using sequence information. Users can chose between four different normalization methods. Please refer to the arguments section below. 
}

\usage{
NormalizeArray(peptideSet, robust = TRUE, standard = FALSE, 
method = "ZpepQuad", centered = TRUE, verbose = FALSE)
}

\arguments{
  \item{peptideSet}{This object contains an ExpressionSet and has the following additional slots}	
  \item{robust}{A \code{logical} value. If TRUE, reweighted least-squares estimates are computed.}
  \item{standard}{A \code{logical}.}
  \item{method}{The normalization method to be used. See Details.}
  \item{centered}{A \code{logical}. If TRUE, recenter the data.}
  \item{verbose}{A \code{logical} value. If verbose is TRUE, progress information is be displayed.}
}

\details{Normalization is applied to each array separately. 

The available methods include \code{"Zpep"}, \code{"ZpepQuad"}, \code{"constant"},
and \code{"binned"}. Methods \code{"Zpep"} and \code{"ZpepQuad"} fit a linear model using either linear or linear and quadratic terms (respectively), regressing
intensity on the peptides' five Z-scale scores. A peptide Z-scale score is obtained by summing over the Z-scale values in Sandburg et al (1998) of the amino acids the peptide comprises.

Peptide Z-scale scores may be provided in the \code{featureRange} slot of \code{peptideSet}. This slot is a \link{RangedData} object \code{x}, and the function will seek five columns labelled \code{z1} through \code{z5} in \code{values(x)}. If these are not found, the function attempts to calculate Z-scales from sequence information found in \code{peptide(peptideSet)}.

Method \code{"constant"} fits a linear model based on amino acid counts of the whole peptide sequence, while method \code{"binned"} uses the amino acid counts of the first, middle, and last thirds of the divided peptide sequence as covariates.

If \code{robust = TRUE} the linear model is fit with \eqn{t_4} distributed errors. If \code{standard = TRUE} probes are grouped into bins based on sorted fitted values, then standardized by the residual standard error within a bin. The method returns the residuals of each peptide intensity in the fitted linear model. If \code{centered = TRUE} the fitted intercept term is added back to the residuals of the fit.
}


\references{
Sandberg, M., Eriksson, L., Jonsson, J., Sjostrom, M., and Wold, S. (1998). New chemical descriptors relevant for the design of biologically active peptides. A multivariate characterization of 87 amino acids. Journal of Medicinal Chemistry 41, 2481-2491.
}

\value{
  A \link{peptideSet} object with updated normalized intensity values.
}

\author{
  Raphael Gottardo, Gregory Imholte
}
