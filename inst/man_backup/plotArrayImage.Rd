\name{plotArrayImage}
\alias{plotArrayImage}

\title{
Plot microarray images
}
\description{
Plot a color image of slide intensities as they are located on the actual
slide. Such plots are helpful for diagnosing spatial abnormalities on slides.
}
\usage{plotArrayImage(peptideSet, array.index = NULL,
  low = "white", high = "steelblue",
  ask = dev.interactive(orNone = TRUE) & 1 < length(array.index))
}
\arguments{
  \item{peptideSet}{
  A \link{peptideSet} object.
}
  \item{array.index}{
  A vector subsetting \code{exprs(peptideSet)}, indicating which slides to plot.
  }
  \item{low}{color of the lowest slide intensity, passed to 
  \code{\link{scale_fill_gradient}}.}

  \item{high}{color of the highest slide intensity, passed to 
  \code{\link{scale_fill_gradient}}.}
  
  \item{ask}{
  logical; if \code{TRUE}, the user is asked before each plot, see
  \code{\link{par}(ask=.)}.
  }
}

\author{
Gregory Imholte
}

\seealso{
\code{\link{plotArrayResiduals}}
}
