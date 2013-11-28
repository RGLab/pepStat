\name{create_db}
\alias{create_db}
\title{Create a peptide collection}
\usage{
  create_db(rd)
}
\arguments{
  \item{rd}{A \code{RangedData} object. The object should
  have a peptide column.}
}
\description{
  Constructor to create peptide collection such as the
  datasets available in PEP.db.
}
\details{
  rd can have additional columns. These columns will be
  kept in the peptide collection.
}
\seealso{
  \code{\link{RangedData}}
}

