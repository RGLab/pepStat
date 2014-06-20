#' Add information to a peptideSet and summarize peptides
#'
#' This function merges the replicates and adds information from a peptide collection
#' to a peptideSet. This collection can include coordinates, alignment information,
#' Z-scales, and other peptide information.
#'
#' @usage summarizePeptides(peptideSet, summary = "median", position = NULL)
#'
#' @param peptideSet A \code{peptideSet}, as created by \code{makePeptideSet}
#' @param summary A \code{character} string. The method used for merging replicates.
#' Available are: "mean" and "median".
#' @param position A \code{data.frame} or \code{GRanges} object. A peptide
#' collection such as the ones available in \code{pepDat}. See details below
#' and vignettes for more information.
#'
#' @return An object of class \code{peptideSet} with added columns and updated ranges.
#'
#' @details
#' The object in the position argument will be passed to \code{create_db}, it
#' can either be a \code{GRanges} object with a peptide as a metadata column, or
#' a \code{data.frame} that can be used to create such \code{GRanges}.
#'
#' Some peptide collections can be found in the \code{pepDat} package.
#'
#'
#' @seealso \code{\link{makePeptideSet}}, \code{\link{create_db}},
#' \code{\link{create_db}}
#'
#' @author Raphael Gottardo, Greory Imholte
#'
#' @rdname summarizePeptides
#'
#' @importFrom GenomicRanges seqnames
#' @export
#' @example examples/pipeline.R
summarizePeptides <- function(peptideSet, summary="median", position=NULL){
	# Check arguments for conformity
	check = .checkArgs_sumPeps(peptideSet, summary, position)
	if(!check){
		stop(attr(check, "ErrorString"))
	}

	df <- as.data.frame(exprs(peptideSet))
	featureSequence <- peptide(peptideSet)

	sdata <- do.call("rbind",
			by(df,list(as.factor(featureSequence)),
					function(x){
						switch(summary,
								mean=colMeans(x, na.rm = TRUE),
								median=rowMedians(t(x), na.rm = TRUE))
					})
	)
	colnames(sdata)<-colnames(df)


	featureID <- sapply(split(featureID(peptideSet),as.factor(featureSequence)),function(x){x[1]})
	featureSequence <- as.character(sapply(split(featureSequence,featureSequence),function(x){x[1]}))

	exprs <- as.matrix(sdata)
	rownames(exprs) <- featureSequence
	colnames(exprs) <- sampleNames(peptideSet)
	nPep <- length(featureID)

	newSet<-new('peptideSet',
			featureRange = GRanges(seqnames = " ", strand = "*",
                             ranges = IRanges(rep(0,nPep),rep(0,nPep)),
					featureID, peptide = featureSequence),
			exprs = as.matrix(sdata),
			experimentData=peptideSet@experimentData)

	sampleNames(newSet) <- sampleNames(peptideSet)


	if(!is.null(position)){
    positiion <- create_db(position)
		# assume that rownames of position GRanges
		# object are peptide sequences in peptideSet,
		# non-null rownames checked in checkArgs above

		# remove elements of GRanges that aren't found in
		# the array
		sub1 <- names(position) %in% peptide(newSet)
		position <- position[sub1,]

		# remove elements of peptideSet that aren't found in
		# GRanges object!
		sub2 <- peptide(newSet) %in% names(position)
		newSet <- newSet[sub2,]

		if(sum(!sub2) > 0){
			message("Some peptides have no match in the GRanges object rownames and are removed from the peptideSet!")
		}

		# reorder peptideSet so that rows of expression matrix
		# match the ordering in the GRanges object
		ind1 <- match(names(position), peptide(newSet))
		newSet <- newSet[ind1,]

    ranges(ranges(newSet)) <- ranges(position)
    values(newSet) <- cbind(values(newSet), values(position))
	}
	pData(newSet) <- pData(peptideSet)
	preproc(newSet)$summary <- summary
	newSet
}

.checkArgs_sumPeps <- function(peptideSet, summary, position){
	OK = TRUE
	attr(OK, "ErrorString") = NULL

	if(!(summary %in% c("median", "mean"))){
		OK = FALSE
		attr(OK, "ErrorString") = ("summary must be either median or mean")
	}
  if(class(peptideSet) != "peptideSet"){
    OK = FALSE
    attr(OK, "ErrorString") = ("peptideSet argument must be an object of class peptideSet")
  }
  return(OK)
}
