\name{clade}
\alias{clade}
\alias{clade-methods}
\alias{clade,RangedData-method}

\title{
	Accessor for the clade column of the pep objects in this package.
}

\description{
	This method returns a \code{"matrix"} of \code{"logical"} with peptides as rownames and clade names as colnames.
}

\section{Methods}{
	\describe{
		\item{\code{signature(object = "RangedData")}}{
			This method is designed for the datasets of this package.}
	}
}

\author{
Renan Sauteraud
}

\examples{
#Get the clades of a dataset
library(PEP.db)
data(pep_hxb2)
clade(pep_hxb2)
}

\keyword{methods}
