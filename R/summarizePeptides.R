summarizePeptides <- function(peptideSet, summary="median", position=NULL)
{    
	if (is(peptideSet,"peptideSet")){
      stop("peptideSet must be an object of class peptideSet")
    }
	
	FUN <- match.fun(summary)
	df <- as.data.frame(exprs(peptideSet))
	featureSequence <- peptide(peptideSet)
	
	sdata<-do.call("rbind", 
	by(df, list(as.factor(featureSequence)), 
	function(x){switch(summary, mean=colMeans(x), median=rowMedians(t(x)))}))
	
	colnames(sdata)<-colnames(df)

	
	featureID <- sapply(split(featureID(peptideSet), as.factor(featureSequence)), function(x){x[1]})
	featureSequence <- as.character(sapply(split(featureSequence, featureSequence), function(x){x[1]}))
		
	exprs <- as.matrix(sdata)
	rownames(exprs) <- featureSequence
	colnames(exprs) <- sampleNames(peptideSet)
	nPep <- length(featureID)

	newSet <- new('peptideSet', 
	featureRange = RangedData(IRanges(rep(0,nPep), 
	rep(0,nPep)), 
	featureID, 
	peptide = featureSequence), 
	exprs = as.matrix(sdata), 
	experimentData = peptideSet@experimentData
	)
		
	sampleNames(newSet) <- sampleNames(peptideSet)

	if (!is.null(position)){
		## sort the positions 
		# Here I assume that the positions are sorted
		# position<-position[order(start(position)),]
		rownames(ranges(newSet)) <- peptide(newSet)
		## Change this line to a match, in case the peptides in the position are not on the array
		newSet <- newSet[rownames(position),]#mainly for subsetting the exprs matrix
		#merge the values from position
		vdata <- c(values(ranges(newSet))[[1]],values(position)[[1]])
		rownames(vdata) <- rownames(values(ranges(newSet)))[[1]]
		values(ranges(newSet)) <- vdata

		#update the IRanges from position
		names(ranges(newSet)) <- names(position)
		ranges(ranges(newSet)) <- ranges(position)
	
		noMatch <- which(!featureSequence%in%rownames(position))
		if (length(noMatch)>0){
			message("Some peptides have no match in the RangedData object and are removed from the peptideSet.")
		}
	}
	pData(newSet) <- pData(peptideSet)
	preproc(newSet)$summary <- summary
	newSet
}
