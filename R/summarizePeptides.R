summarizePeptides<-function(peptideSet,summary="mean",position=NULL,...)
{

  FUN <- match.fun(summary)
  df<-as.data.frame(exprs(peptideSet))
  featureSequence<-peptide(peptideSet)
	
#	system.time(
			sdata<-do.call("rbind"
							,by(df,list(as.factor(featureSequence))
									,function(x){
#										browser()
										switch(summary ,
												mean=colMeans(x),
												median=rowMedians(t(x))
												)			
										
									})
							)
#	)
			colnames(sdata)<-colnames(df)

	
#	featureSequence<-sapply(split(peptideSet@featureSequence,as.factor(peptideSet@featureSequence)),unique)
#	featureID<-sapply(split(peptideSet@featureID,as.factor(peptideSet@featureSequence)),unique)

# Since they are the same, we only take the first one
	featureID<-sapply(split(featureID(peptideSet),as.factor(featureSequence)),function(x){x[1]})
	featureSequence<-as.character(sapply(split(featureSequence,featureSequence),function(x){x[1]}))
		
	exprs<-as.matrix(sdata)
	rownames(exprs)<-featureSequence
	colnames(exprs)<-sampleNames(peptideSet)
	nPep<-length(featureID)
	newSet<-new('peptideSet'
				,featureRange=RangedData(IRanges(rep(0,nPep),rep(0,nPep))
										,featureID								
										,peptide=featureSequence
										)
				, exprs=as.matrix(sdata)
				, experimentData=peptideSet@experimentData
			)
		
	sampleNames(newSet)<-sampleNames(peptideSet)

#	browser()  
	
  if(!is.null(position))
  {
	##sort the positions    
	position<-position[order(start(position)),]
	rownames(ranges(newSet))<-peptide(newSet)
	newSet<-newSet[rownames(position),]#mainly for subsetting the exprs matrix
#	browser()
	#merge the values from position
	vdata<-c(values(ranges(newSet))[[1]],values(position)[[1]])
	rownames(vdata)<-rownames(values(ranges(newSet)))[[1]]
	values(ranges(newSet))<-vdata

	#update the IRanges from position
	names(ranges(newSet))<-names(position)
	ranges(ranges(newSet))<-ranges(position)
	
	noMatch<-which(!featureSequence%in%rownames(position))
	if(length(noMatch)>0)
	{
		message(featureSequence[noMatch]," have no position provided and are removed from the peptideSet!")
	}

  }
  
#  if(!is.null(annotation))
#  {
#    annotation<-annotation[!duplicated(as.character(annotation$Sequence)),]
#    match.seq<-match(newSet@featureSequence,annotation$Sequence)
#    annotation<-annotation[match.seq,]
#    newSet@featureAnnotation<-annotation
#  }

#	newSet[is.finite(newSet@featurePosition),]
	newSet
}
