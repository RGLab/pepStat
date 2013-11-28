\name{plotArrayImage}
\alias{plotArrayImage}
\alias{plotArrayResiduals}
\title{Plot microarray images}
\usage{
  plotArrayImage(peptideSet, array.index = NULL,
    low = "white", high = "steelblue",
    ask = dev.interactive(orNone = TRUE) & 1 < length(array.index))

  plotArrayResiduals(peptideSet, array.index = NULL,
    smooth = FALSE, low = "blue", high = "red",
    ask = dev.interactive(orNone = TRUE) & 1 < length(array.index))
}
\arguments{
  \item{peptideSet}{A \code{peptideSet} object}

  \item{array.index}{A vector subsetting
  \code{exprs(peptideSet)}, indicating which slides to
  plot}

  \item{smooth}{A \code{logical}, a 2D spatial smoother is
  applied to residuals, the fitted residuals are plotted.}

  \item{low}{A \code{character} string. The color of the
  lowest slide intensity. passed to
  \code{scale_fill_gradient}. the fitted residuals are
  plotted.}

  \item{high}{A \code{character} string. The color of the
  highest slide intensity. passed to
  \code{scale_fill_gradient}.}

  \item{ask}{A \code{logical}. If TRUE, the user is asked
  before each plot. See \code{par(ask=.)}.}
}
\description{
  Plot a color image of the intensities on a slide. These
  plots are helpful to diagnose spatial abnormalities.
}
\author{
  Gregory Imholte
}

