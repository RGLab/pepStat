#now we can only assume there is one chromosome per time
#Furthermore, the vectors must be sorted

slidingMean<-function(peptideSet, width=5, verbose=FALSE)
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
  
  #whether the Control array is available
  y<-exprs(peptideSet)
  
  if(length(peptideSet)==0)
  {
    stop("A position vector must be specified in the peptideSet")
  }

  # This could be made more efficient with multicore
  ny<-apply(y,2,slidingMeanVector,width,position(peptideSet))

  # y<-matrix(obj$SScore,nPep,ncol(y),byrow=TRUE)
  exprs(peptideSet)<-ny
  
  if(verbose)
  {
    cat("** Finished processing ",nrow(y)," probes on ",ncol(y)," arrays **\n");
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