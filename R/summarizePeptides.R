summarizePeptides <- function(peptideSet, summary="median", position=NULL)
{
	# Check arguments for conformity
	check = .checkArgs_sumPeps(peptideSet, summary, position)
	if(!check)
		stop(attr(check, "ErrorString"))
	
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
			featureRange = RangedData(IRanges(rep(0,nPep),rep(0,nPep)),
					featureID, peptide = featureSequence),
			exprs = as.matrix(sdata),
			experimentData=peptideSet@experimentData)
	
	sampleNames(newSet) <- sampleNames(peptideSet)


	if(!is.null(position))
	{
		# assume that rownames of position RangedData 
		# object are peptide sequences in peptideSet,
		# non-null rownames checked in checkArgs above 

		# remove elements of RangedData that aren't found in 
		# the array
		sub1 <- rownames(position) %in% peptide(newSet)
		position <- position[sub1,]
		
		# remove elements of peptideSet that aren't found in
		# RangedData object!
		sub2 <- peptide(newSet) %in% rownames(position)
		newSet <- newSet[sub2,]
		
		if(sum(!sub2) > 0){
			message("Some peptides have no match in the RangedData object rownames and are removed from the peptideSet!")
		}
		
		# reorder peptideSet so that rows of expression matrix
		# match the ordering in the RangedData object
		ind1 <- match(rownames(position), peptide(newSet))
		newSet <- newSet[ind1,]
		
		#split and reorder peptide info from newSet so that
		# order of spaces matches position
		DframeTmp <- split(values(ranges(newSet))[[1]], space(position))
		ind2 <- match(names(values(position)), names(DframeTmp))
		
		DframeTmp <- DframeTmp[ind2]
		rownames(DframeTmp) <- rownames(values(position))
		
		if(!all(names(DframeTmp) == names(values(position))))
			stop("mismatched RangedData objects")
		
		vdata <- cbind(DframeTmp, values(position))
		rdata <- ranges(position)
		
		newSet@featureRange <- RangedData(rdata, vdata)
	}
	pData(newSet) <- pData(peptideSet)
	preproc(newSet)$summary <- summary
	newSet
}

.checkArgs_sumPeps = function(peptideSet, summary, position)
{
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
	
	if(class(position) != "RangedData" & !is.null(position)){
		OK = FALSE
		attr(OK, "ErrorString") = ("if non-NULL, position argument must be a RangedData object")
	}
	
	if(is.null(rownames(position)) & !is.null(position)){
		OK = FALSE
		attr(OK, "ErrorString") = ("rownames of position RangedData object must be peptide sequences")
	}
	
	return(OK)	
}
