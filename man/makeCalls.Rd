\name{makeCalls}
\alias{makeCalls}
\title{Make antibody binding positivity calls}
\usage{
makeCalls(peptideSet, cutoff = 1.2, method = "absolute", freq = TRUE,
  group = NULL, verbose = FALSE)
}
\arguments{
  \item{peptideSet}{A \code{peptideSet} object. The
  peptides, after normalization and possibly data
  smoothing.}

  \item{cutoff}{A \code{numeric}. If FDR, the FDR
  threshold. Otherwise, a cutoff for the background
  corrected intensities.}

  \item{method}{A \code{character}. The method used to make
  positivity calls. "absolute" and "FDR" are available. See
  details below.}

  \item{freq}{A \code{logical}. If set to TRUE, return the
  percentage of slides calling a peptide positive.
  Otherwise, return a \code{logical} indicating binding
  events.}

  \item{group}{A \code{character}. Only used when freq is
  set to TRUE. A character indicating a variable by which
  to group slides. If non-null the percentage is calculated
  by group.}

  \item{verbose}{A \code{logical}. If set to TRUE, progress
  information will be displayed.}
}
\value{
If freq = TRUE, a \code{numeric} \code{matrix} with
peptides as rows and groups as columns where the values are
the frequency of response in the group. If freq = FALSE, a
\code{logical} \code{matrix} indicating binding events for
each peptide in each subject.
}
\description{
After normalization and data smoothing, this last step
makes the call for each peptide of the peptideSet after
baseline correcting the peptide intenstities.
}
\details{
This function requires specific variables ptid and visit in
pData(peptideSet). The variable \code{ptid} should indicate
subjects, and the variable \code{visit} should be a factor
with levels pre and post.

If slides are paired for subjects, intensities
corresponding to post-visit are substracted from pre. If
slides are not paired, slides with pre have intensities
averaged by peptides, and averaged peptide intensities are
subtracted from slides that have entry post. Calls are made
on these baseline corrected intensities.

When method = FDR, a left-tail method is used to generate a
threshold controlling the False Discovery Rate at level
\code{cutoff}. When method = absolute, Intensities
exceeding the threshold are labelled as positive.

When freq = TRUE a group variable may be specified. The
argument group indicates the name of a variable in
pData(peptideSet) by which positive calls should be
grouped. The call frequency for each peptide is calculated
within groups.
}
\author{
Greg Imholte
}

