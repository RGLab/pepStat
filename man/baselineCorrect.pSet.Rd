\name{baselineCorrect.pSet}
\alias{baselineCorrect.pSet}
\title{Substract baseline intensities}
\usage{
  baselineCorrect.pSet(pSet, verbose = FALSE)
}
\arguments{
  \item{pSet}{A \code{peptideSet} with sample PRE and POST
  visits.}
}
\description{
  Correct intensities by substracting PRE visit sample
  intensities.
}
\details{
  If samples are not PAIRED (One PRE and POST for each
  ptid), then the average expression of all PRE visit
  samples is substracted from each sample.
}

