\name{plotArrayResiduals}
\alias{plotArrayResiduals}

\title{
Plot image of residuals for microarrays
}
\description{
Plot a color image of residual probe intensities as they are located on the
actual slide. Such plots are helpful for diagnosing spatial abnormalities on
slides.
}
\usage{plotArrayResiduals(peptideSet, array.index = NULL, smooth = FALSE,
  low = "blue", high = "red",
  ask = dev.interactive(orNone = TRUE) & 1 < length(array.index))
}
\arguments{
  \item{peptideSet}{
  A \link{peptideSet} object.
}
  \item{array.index}{
  A vector subsetting \code{exprs(peptideSet)}, indicating which slides to plot.
  }
  \item{smooth}{logical; if \code{TRUE}, a 2d spatial smoother is applied to 
  residuals, and the fitted residuals are plotted.}
  \item{low}{color of the lowest slide intensity, passed to 
  \code{\link{scale_fill_gradient}}.}

  \item{high}{color of the highest slide intensity, passed to 
  \code{\link{scale_fill_gradient}}.}
  
  \item{ask}{
  logical; if \code{TRUE}, the user is asked before each plot, see
  \code{\link{par}(ask=.)}.
  }
}
\details{
 Residuals are based on within-slide replications of the same type of probe.
 A probe's residual is its original intensity 
 minus the mean of intensities among probes of the same type.
 
 When \code{smooth} is \code{TRUE}, the \code{\link{smooth.2d}} function in 
 the \code{fields} package is applied to the residuals, possibly aiding in the
 detection of spatial effects.
}

\author{
Gregory Imholte
}

\seealso{
\code{\link{plotArrayImage}}
}