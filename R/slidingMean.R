#now we can only assume there is one chromosome per time
#Furthermore, the vectors must be sorted

slidingMean <-function(peptideSet, width=5, verbose=FALSE, split.by.space=TRUE)
{
	if(class(peptideSet)!="peptideSet")
	{
		stop("peptideSet must be an object of class peptideSet")
	}
	if(preproc(peptideSet@experimentData)$transformation!="log" & preproc(peptideSet@experimentData)$transformation!="vsn")
	{
		stop("The probe measurements need to be log/vsn transformed!")
	}
	if(preproc(peptideSet@experimentData)$normalization=="none")
	{
		warning("You should probably normalize your data before using this function")
	}
	
	# Check whether user wishes to apply to multiple spaces 
	if(split.by.space & length(names(ranges(peptideSet))) > 1)
	{
    oldrownames = rownames(ranges(peptideSet))
		s <- space(ranges(peptideSet))
		pSetList <- split(peptideSet, s)
		smoothedpSetList <- lapply(pSetList, function(set)
        {
          p <- position(set)
					o <- order(p)
					
					# reorder peptideSet by position order
					set <- set[o,]
				
					y <- exprs(set)
					p <- positon(set)
					ny <- apply(y, 2, slidingMeanVector, width, p)
					exprs(set) <- ny
					rownames(exprs(set)) <- rownames(y)
					
					set
				})
		
		# collect the sorted RangedData into a single RangedData object,
		# put back into peptideSet
		rdata = do.call("rbind", lapply(smoothedpSetList, ranges))
		peptideSet@featureRange = rdata
		
		# collect the separate smoothed intensities into a single expression
		# matrix, gather rownames
		exprs(peptideSet) = do.call("rbind", lapply(smoothedpSetList, exprs))
		rownames(exprs(peptideSet)) = do.call("c", 
				lapply(smoothedpSetList, function(x) rownames(exprs(x))))
	}
	
  else
	{
		# This could be made more efficient with multicore
		p <- position(peptideSet)
		o <- order(p)
		ro <- order(o)
		
		y <- exprs(peptideSet)[o,]
    p <- position(peptideSet)[o]
		ny <- apply(y, 2, slidingMeanVector, width, p)
		exprs(peptideSet) <- ny[ro,]
	}
	
	if(verbose)
	{
		cat("** Finished processing ",nrow(peptideSet)," probes on ",ncol(peptideSet)," arrays **\n");
	}
	peptideSet
}

slidingMeanVector<-function(y,dMax,position)
{
	obj<-.C("MATScore",
			y,
			as.integer(length(y)),
			as.integer(1),
			as.integer(position),
			as.double(dMax),
			SScore=double(length(y)),
			as.integer(rep(1,length(y))),
			package="pepStat")
	
	return(obj$SScore);
}

callEnrichedRegions<-function(SScore, dMax=2, dMerge=2, nProbesMin=2, method="score", threshold=5, verbose=FALSE)
{
	if(!is(SScore,"RangedData"))
	{
		stop("SScore must be an object of class RangedData (e.g. returned by computeSScore)!")
	}
	
	if(method=="score")
	{
		methodNum<-1
	}
	else if(method=="pValue")
	{
		if(threshold>1 | threshold<0)
		{
			stop("When the pValue method is selected, the threshold should be between 0 and 1")
		}
		methodNum<-2
	}
	else if(method=="FDR")
	{
		if(threshold >1 | threshold<0)
		{
			stop("When the FDR method is selected, the threshold should be between 0 and 1")
		}
		methodNum<-3
	}
	else
	{
		stop("Argument 'Method' must be either score or pValue or FDR")
	}
	chrAll<-space(SScore)
	seqNum<-as.numeric(as.factor(chrAll))
	nProbes<-length(seqNum)
	startMat<-start(SScore)
	scoreMat<-SScore$score
	
	obj<-.C("callEnrichedRegions",
			as.double(scoreMat),
			as.integer(nProbes),
			as.integer(startMat),
			as.double(dMerge),
			as.double(dMax),  
			as.double(threshold),
			pValue=double(nProbes),
			as.integer(methodNum),
			regions=integer(nProbes), 
			as.integer(verbose),
			as.integer(seqNum),
			numRegions = integer(1),
			package="pepStat")
	
	pValue<-obj$pValue
	
	if(verbose)
	{
		cat("** Number of Enriched regions is ", obj$numRegions, " **\n")
	}
	numRegions<-obj$numRegions
	if(numRegions > 0)
	{
		# We have detected some regions
		obj<-.C("getIndices",
				as.integer(obj$regions),
				as.integer(nProbes),
				as.integer(numRegions),
				Start=integer(numRegions),
				End=integer(numRegions),
				package="pepStat")
		
		center<-score<-start<-end<-chr<-rep(0,numRegions)
		
		for(i in 1:numRegions)
		{
			start[i]<-startMat[obj$Start[i]]
			end[i]<-startMat[obj$End[i]]
			score[i]<-max(scoreMat[obj$Start[i]:obj$End[i]])
			center[i]<-startMat[(obj$Start[i]:obj$End[i])[which.max(scoreMat[obj$Start[i]:obj$End[i]])]]
			chr[i]<-chrAll[obj$Start[i]]
		}
		
		ind<-obj$End-obj$Start+1>=nProbesMin
		
		if(verbose)
		{
			cat("** ", sum(!ind), " enriched regions have less than ", nProbesMin, " probes and are filtered out**\n")
		}
		ranges<-IRanges(start=start[ind],end=end[ind])
		chr<-chr[ind]
		score<-score[ind]
		center<-center[ind]
	}
	else
	{
		if(verbose)
		{
			cat("** No regions to output **\n")
		}
		ranges<-IRanges(NULL)
		return(RangedData(ranges))
	}
	# Here I assume that all the sequences have the same length
	RD<-RangedData(ranges, space = chr, score=score, center=center)
	RD
}
