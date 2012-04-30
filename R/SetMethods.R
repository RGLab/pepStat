setMethod("show", "peptideSet",function(object){
    cat("Object of class 'peptideSet' contains","\n")
	print(as(object,"ExpressionSet"))
	print(ranges(object))
#    cat("featureSequence, featureAnnotation\n")
})

setAs(from="peptideSet",to="ExpressionSet",function(from){
#			browser()
			ExpressionSet(assayData(from)
					,phenoData=phenoData(from)
#					,featureData=featureData(from)
					,experimentData=experimentData(from)
					,annotation=annotation(from)
					,protocolData=protocolData(from)
					)
})

setMethod("summary", signature("peptideSet"),
    function(object) {
      cat("   Sample name(s): ",sampleNames(object@phenoData)," \n")
      cat("   The total number of probes is: ",length(peptide(object))," \n")
      cat("   Preprocessing Information \n")
      cat("     - Transformation:",preproc(object@experimentData)$transformation, "\n")
      cat("     - Normalization:",preproc(object@experimentData)$normalization, "\n")      
    }
)

setMethod("[","peptideSet",
		function(x,i,j,..., drop=FALSE)
		{
			if(!missing(i))
			{
			sdata<-exprs(x)[i,j]		
			featureRange=ranges(x)[i,]
			}
			else
			{
			    sdata<-exprs(x)[,j]				
				featureRange<-ranges(x)
			}
		    newSet<-new('peptideSet'
							,exprs=as.matrix(sdata)
							,featureRange=featureRange
							,experimentData=x@experimentData)
			pData(newSet)<-pData(x)[j,]
			sampleNames(newSet)<-sampleNames(x)[j]
		    newSet
			
		})


setGeneric("position", function(x, ...) standardGeneric("position"))
setMethod("position","peptideSet",function(x){
			#x@featurePosition
			round((start(ranges(x))+end(ranges(x)))/2)
		})


setMethod("start","peptideSet",function(x){
			#x@featurePosition
			start(ranges(x))
		})

setMethod("end","peptideSet",function(x){
			#x@featurePosition
			end(ranges(x))
		})


setReplaceMethod("ranges","peptideSet",
		function(x,value)
		{
			x@featureRange<-value
			x
		})

setMethod("ranges","peptideSet",
		function(x)
		{
			x@featureRange
		})

setGeneric("peptide", function(x, ...) standardGeneric("peptide"))
setMethod("peptide","peptideSet",
		function(x,type=NULL)
		{
			#x@featureSequence
			validTypes<-c("peptide","aligned","trimmed")
			if(is.null(type))
			{
				type<-"peptide"
				
			}
			
			
			if(type%in%validTypes)
			{
				ranges(x)[[type]]	
			}else
			{
				warning("'",type, "' is not valid types(",validTypes,")!")
			}
		})

setGeneric("featureID", function(x, ...) standardGeneric("featureID"))
setMethod("featureID","peptideSet",
		function(x,type=NULL)
		{
			#x@featureSequence
			if(is.null(type))
			{
				ranges(x)[["featureID"]]
			}
		})


setGeneric("clade", function(x, ...) standardGeneric("clade"))

setMethod("clade","peptideSet",function(x){
#			browser()
#	x@featureAnnotation[names(x@featureAnnotation)%in%c("M","A","B","C","D","CRF01","CRF02")]
	df<-IRanges::as.data.frame(values(ranges(x)))
	df[,colnames(df)%in%c("M","A","B","C","D","CRF01","CRF02")]		
	
})

setGeneric("split", function(x,f, ...) standardGeneric("split"))

setMethod("split","peptideSet",function(x, f, byrow=TRUE){
  if(is.vector(f) | is.factor(f))
  {
    f<-as.factor(f)
    if(byrow)
    {
      lapply(1:nlevels(f),function(i,pSet,f){pSet[f==levels(f)[i],]},x,f)
    }
    else
    {
      lapply(1:nlevels(f),function(i,pSet,f){pSet[,f==levels(f)[i]]},x,f)
    }
  }
  else
  {
    lapply(1:ncol(f),function(i,pSet,f){pSet[f[,i],]},x,f)
  }
})


## Bind several tiling sets together
#
#setMethod("rbind", "tilingSet", function(..., deparse.level=1) {
#  args <- list(...)
#
#  names<-lapply(args,function(x){sampleNames(x@phenoData)})
#  lnames<-unlist(lapply(names,"length"))
#  # Check that the number of names is the same for all
#  if(length(unique(lnames))!=1)
#  {
#    stop("Objects to be concatenated should have the same number of columns")    
#  }
#
#  featureChromosome<-unlist(lapply(args,function(x){x@featureChromosome}))
#  featurePosition<-unlist(lapply(args,function(x){x@featurePosition}))
#  featureCopyNumber<-unlist(lapply(args,function(x){x@featureCopyNumber}))
#  featureSequence<-unlist(lapply(args,function(x){x@featureSequence}))
#  
#  ord<-order(featureChromosome,featurePosition)
#  
#  featureChromosome<-unlist(lapply(args,function(x){x@featureChromosome}))
#  y<-do.call("rbind",lapply(args,function(x){exprs(x)}))
#  
#  ## For the sample name, I simply use the first set
#  colnames(y)<-names[[1]]
#  rownames(y)<-NULL
#  newSet<-new('tilingSet', featureChromosome=featureChromosome[ord],featurePosition=featurePosition[ord],
#  featureCopyNumber=featureCopyNumber[ord], exprs=y[ord,], genomeName=unlist(lapply(args,function(x){x@genomeName})),
#  featureSequence=featureSequence[ord], experimentData=args[[1]]@experimentData)
#  return(newSet)
#})
#
