\name{makePeptideSet}
\alias{makePeptideSet}

\title{peptideSet constructor}

\description{
This function reads GenePix results (.gpr) files and creates a peptideSet object gathering experiment information.
}

\usage{
makePeptideSet(files = NULL, path = NULL, mapping.file = NULL, 
  use.flags = FALSE, rm.control.list = NULL, norm.empty = FALSE,
	empty.control.list = c("empty","blank control"),
  bgCorrect.method = "normexp", log = TRUE, check.row.order = FALSE,
  verbose = FALSE)
}

\arguments{
\item{files}{A character vector containing the names of the files to read. If omitted, all the files with the extension .gpr in specified path will be read.}
\item{path}{A character string giving the directory containing files.}
\item{mapping.file}{A character string or a data frame giving a mapping file for the data. See Details.}
\item{use.flags}{A \code{logical}. Should spots with spot flag values -99 or lower be excluded?}
\item{rm.control.list}{A vector of character strings. The name of the controls to be excluded from the peptideSet.}
\item{norm.empty}{A \code{logical}. Subtract mean intensity of empty spots from remaining intensities within slide?}
\item{empty.control.list}{A vector of character strings. The name of the empty controls. These are used or background correction if \code{norm.empty = TRUE}.}
\item{bgCorrect.method}{A character string. The name of the method used for background correction. See Details.}
\item{log}{A \code{logical}. Should intensities be \code{log2} transformed after BG correction?}
\item{check.row.order}{A \code{logical}. Should slides be reduced to a common set of peptides?}
\item{verbose}{A \code{logical}.}
}

\details{
	GenePix results files (\code{.gpr}) are read when found in either the \code{files} or \code{path} arguments. By default, the foreground and background median intensities of the "red" channels, "F635 Median" and "B635 Median", are read. The background correction specified in \code{bgCorrect.method} is passed to the \link{backgroundCorrect} method in the \code{limma} package.

	The argument \code{mapping.file} may either be a character string file path or data frame. If a file path, it must point to a \code{.csv} file which is read and used as experiment information. Each row of the \code{mapping.file} data frame is taken to correspond with information about a single slide. Entries in the data frame are matched to arrays by file name. One column in the \code{mapping.file} data frame MUST be labeled "Filename" (not case sensitive) if provided. Entries in this column must match file names read from the \code{path} or \code{files} arguments after dropping the \code{.gpr} extension. For example, "PID_101_vst1.gpr" would be matched to "PID_101_vst1" in a "filename" column provided by \code{mapping.file}.
	
	The argument \code{rm.control.list} specifies (possibly partial) names of array probes to be removed from the data set. Elements of this vector are matched against the Name and ID columns of the \code{.gpr} files.
	
	After background correction, if \code{norm.empty = TRUE}, intensities of probes matching strings (possibly partially) in \code{empty.control.list} are averaged within each array, and the average is subtracted from all intensities within the each array.
	
	If \code{check.row.order = TRUE}, the final set of probes is taken to be those with IDs found in all arrays that were read.
	
	Row, Column, and Block spatial array positions for each probe are stored in the \code{featureRanges} slot of the returned object.
}

\value{
	A \link{peptideSet} object.
}

\author{Raphael Gottardo, Gregory Imholte}

\seealso{
\link{peptideSet}
}

\keyword{functions}
