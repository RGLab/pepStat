\name{makeCalls}
\alias{makeCalls}

\title{Make antibody binding positivity calls}

\description{
After normalization and data smoothing, this last step makes the calls for each peptide of the peptideSet after baseline correcting the peptide intensities.
}

\usage{
makeCalls(peptideSet, cutoff = .1, method = "absolute", freq = TRUE, group = NULL, verbose = FALSE)
}

\arguments{
  \item{peptideSet}{A \code{peptideSet} object. The peptides, after normalization and possibly data smoothing.}
  \item{cutoff}{A \code{numeric}.}
  \item{method}{A \code{character}. The method used to make the positivity calls. See Details.}
	\describe{
		\item{FDR}{Threshold intensities based on controlling the False Discovery Rate (FDR). \code{cutoff} indicates the desired FDR.}
		\item{absolute}{A simple threshold. Threshold values using \code{cutoff}.}
		}
  \item{freq}{A \code{logical}. Return percentage of slides calling a peptide positive? If \code{group} is non-null, the percentage is calculated by group. See details.}
  \item{group}{A \code{character}. Only used when \code{freq = TRUE}. A character indicating a variable by which to group slides.}
  \item{verbose}{A \code{logical}. If set to TRUE, progress information will be displayed.}
}

\details{
This function requires specific variables \code{ptid} and \code{visit} in \code{pData(peptideSet)}. The variable \code{ptid} should indicate subjects, and the variable \code{visit} should be a factor with levels \code{pre} and \code{post}. 

If slides are paired for subjects, intensities corresponding to \code{post} are subtracted from \code{pre}. If slides are not paired, slides with \code{pre} have intensities averaged by peptides, and averaged peptide intensities are subtracted from slides that have entry \code{post}. Calls are made on these baseline corrected intensities.

When \code{method = FDR}, a left-tail method is used to generate a threshold controlling the False Discovery Rate at level \code{cutoff}. Intensities exceeding the threshold are labelled as positive.

When \code{freq = TRUE} a \code{group} variable may be specified. The argument \code{group} indicates the name of a variable in \code{pData(peptideSet)} by which  positive calls should be grouped. The call frequency for each peptide is calculated within groups. 
}

\value{
If \code{freq = TRUE}, a \code{numeric} vector with one value for each peptide and each \code{group}, if specified. The percentage of calls made within each group for each peptide.

If \code{freq = FALSE}, a \code{logical} matrix with TRUE values indicating a positive call for a subject. Rows correspond to peptides, columns correspond to subjects in \code{pData(peptideSet)$ptid}.
}


\seealso{
\code{NormalizeArray}
\code{slidingMean}
}

\keyword{functions}
