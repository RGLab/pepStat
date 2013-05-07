setMethod("show", "peptideSet",function(object){
  cat("Object of class 'peptideSet' contains","\n")
  print(as(object,"ExpressionSet"))
  print(ranges(object))
})

setAs(from="peptideSet",to="ExpressionSet",function(from){
			ExpressionSet(assayData(from)
					,phenoData=phenoData(from)
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

setMethod("[", signature("peptideSet", i = "ANY", j = "ANY"),
          function(x, i, j, ..., drop = FALSE) {
            if (!missing(i)) {
              sdata <- exprs(x)[i, j, drop = drop]
              featureRange <- ranges(x)[i, ]
            } 
            else {
              sdata <- exprs(x)[, j, drop = drop]
              featureRange <- ranges(x)
            }
            
            newSet<-new('peptideSet', 
                        exprs = as.matrix(sdata), 
                        featureRange = featureRange, 
                        experimentData = x@experimentData)
            
            if (!missing(j)) {
             pData(newSet) <- pData(x)[j,]
             sampleNames(newSet) <- sampleNames(x)[j]
            } else {
              pData(newSet) <- pData(x)
            }  
            newSet
          })

setMethod("subset", signature(x = "peptideSet"), 
          function (x, subset, drop = FALSE, ...) {
            if (missing(subset)){
              r <- rep(TRUE,nrow(x)) 
            } 
            else {
              e <- substitute(subset)#class(e) = call
              r <- eval(e, pData(x), parent.frame())
              if (!is.logical(r)) 
                stop("'subset' must evaluate to logical")
              r <- r & !is.na(r)
              vars<-r
            }
            x[, vars, drop = drop]
          })


setGeneric("position", function(x, ...) standardGeneric("position"))
setMethod("position", "peptideSet", function(x){
			round((start(ranges(x))+end(ranges(x)))/2)
		})


setMethod("start", "peptideSet", function(x){
			start(ranges(x))
		})

setMethod("end", "peptideSet", function(x){
			end(ranges(x))
		})

setMethod("width", "peptideSet", function(x){
			width(ranges(x))
		})

setMethod("values", "peptideSet", function(x){
			values(ranges(x))
		})

setMethod("values<-", "peptideSet", function(x, value){
			values(ranges(x)) <- value
			return(x)
		})

setReplaceMethod("ranges", "peptideSet",
		function(x,value)
		{
			x@featureRange <- value
			x
		})

setMethod("ranges", "peptideSet",
		function(x)
		{
			x@featureRange
		})

setGeneric("peptide", function(x, ...) standardGeneric("peptide"))
setMethod("peptide", "peptideSet",
		function(x, type=NULL)
		{
			validTypes<-c("peptide", "aligned", "trimmed")
			if (is.null(type)){
				type <- "peptide"
			}
			
			if (type%in%validTypes){
				ranges(x)[[type]]	
			} 
      else {
				warning("'",type, "' is not valid types(",validTypes,")!")
			}
		})

setGeneric("featureID", function(x, ...) standardGeneric("featureID"))
setMethod("featureID", "peptideSet",
		function(x, type=NULL){
			if (is.null(type))
			{
				ranges(x)[["featureID"]]
			}
		})

setGeneric("split")
setMethod("split", "peptideSet", function(x, f, byrow = TRUE){
  if(is.vector(f) | is.factor(f))
  {
    f <- as.factor(f)
    if(byrow)
    {
      lapply(1:nlevels(f), 
             function(i, pSet, f){pSet[f == levels(f)[i],]}, x, f)
    }
    else
    {
      lapply(1:nlevels(f),
             function(i, pSet, f){pSet[, f == levels(f)[i]]}, x, f)
    }
  }
  else
  {
    lapply(1:ncol(f), function(i, pSet, f){pSet[f[, i], ]}, x, f)
  }
})


setMethod("cbind", "peptideSet", function(..., deparse.level=1){
 args <- list(...)

 names <- unlist(sapply(args, function(x){sampleNames(x)}))
 pd.list <- lapply(args, function(x){pData(x)})
 pd<-do.call(rbind, pd.list)
 
 eSet.list <- lapply(args, function(x){exprs(x)})
 eSet <- do.call(cbind, eSet.list)
 
 newSet<-new('peptideSet',
             exprs = as.matrix(eSet),
             featureRange = args[[1]]@featureRange,
             experimentData = args[[1]]@experimentData)
 pData(newSet) <- pd
 sampleNames(newSet) <- names
 newSet 
})

setGeneric("write.pSet", 
           function(x, bg.correct=FALSE, ...) standardGeneric("write.pSet"))
setMethod("write.pSet", "peptideSet", function(x, bg.correct=FALSE, ...){
  if (bg.correct) {
    exprs<-pepStat:::.bgCorrect.pSet(x)
  } else {
    exprs<-exprs(x)
  }
  y <- cbind(peptide(x), start(x), end(x), featureID(x), exprs)
  colnames(y)[1:4] <- c("peptide", "start", "end", "annotation")
  write.csv(y, ...)
})
	

#clade acessor
setGeneric("clade",
                def = function(object)
                        standardGeneric("clade"))

setMethod("clade",
          signature = signature(object = "RangedData"),
          definition = function(object){
            cladeList <- unique(unlist(
              strsplit(levels(as.factor(object$clade)), ","))) #List of all possible clades
            len <- nrow(object)
            retMatrix <- matrix(FALSE, nrow = len, ncol = length(cladeList))
            pepClades <- strsplit(object$clade, split = ",") #clades for each peptide
            for(pepIdx in 1:len){
              tmpList <- cladeList %in% pepClades[[pepIdx]]
              retMatrix[pepIdx, ] <- tmpList
            }   
            rownames(retMatrix) <- rownames(object)
            colnames(retMatrix) <- cladeList
            return(retMatrix)
          })  

setMethod("clade",
                signature = signature(object="peptideSet"),
                definition = function(object){ clade(object@featureRange)
                }) 
