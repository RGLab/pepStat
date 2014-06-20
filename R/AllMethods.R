#' peptideSet methods
#'
#' Methods for handling peptideSet objects
#' @name peptideSet-methods
#' @rdname peptideSet-methods
#'
#' @section Accessors:
#' \describe{
#'  \item{\code{nrow(x)}:}{The number of peptides in x.}
#'  \item{\code{ncol(x)}:}{The number of samples in x.}
#'  \item{\code{start(x)}:}{Get the starts of the peptides.}
#'  \item{\code{end(x)}:}{Get the ends of the peptides.}
#'  \item{\code{width(x)}:}{Get the widths of the peptides.}
#'  \item{\code{position(x)}:}{Get the coordinates of the central amino-acid of
#'  each peptide, given by: \code{round((start(x) + end(x))/2)}.}
#'  \item{\code{ranges(x)}:}{Returns a \code{GRanges} object that contains
#'  the annotations for the peptides.}
#'  \item{\code{ranges(x)<- value}}{Set annotations for the peptides.}
#'  \item{\code{values(x)}:}{Returns a \code{SplitDataFrameList}. Accessor for the
#'  values of the featureRange slot.}
#'  \item{\code{clade(x)}:}{If available, returns the clade information for each
#'  peptide as a \code{matrix}.}
#'  \item{\code{peptide(x)}:}{Get the sequence of the peptides.}
#'  \item{\code{peptide(x) <- value}}{Set the sequence of the peptides.}
#'  \item{\code{featureID(x)}:}{Get the ID of the peptides.}
#'  \item{\code{pepZscore(x)}:}{If available, returns a \code{matrix} of the zScores
#'  for each peptide.}
#'  \item{\code{pepZscore(x) <- value}}{Set the zScores for each peptide}
#' }
#'
#' @aliases
#' start,peptideSet-method
#' end,peptideSet-method
#' width,peptideSet-method
#' position
#' position-method
#' position,peptideSet-method
#' ranges,peptideSet-method
#' ranges<-,peptideSet-method
#' values,peptideSet-method
#' clade
#' clade-methods
#' clade,GRanges-method
#' clade,peptideSet-method
#' peptide
#' peptide<-
#' peptide-method
#' peptide,peptideSet-method
#' peptide<-,peptideSet,character-method
#' featureID
#' featureID-method
#' featureID,peptideSet-method
#' pepZscore
#' pepZscore<-
#' pepZscore-method
#' pepZscore,peptideSet-method
#' pepZscore<-,peptideSet,data.frame-method
#' pepZscore,GRanges-method
#' pepZscore<-,GRanges,data.frame-method
#' [,peptideSet,ANY,ANY,ANY-method
#' subset,peptideSet-method
#' show,peptideSet-method
#' summary,peptideSet-method
#'
#' @exportMethod "start"
#' @exportMethod "end"
#' @exportMethod "width"
#' @exportMethod "position"
#' @exportMethod "ranges"
#' @exportMethod "values"
#' @exportMethod "ranges<-"
#' @exportMethod "clade"
#' @exportMethod "peptide"
#' @exportMethod "pepZscore"
#' @exportMethod "featureID"
#' @exportMethod "peptide<-"
#' @exportMethod "pepZscore<-"
#'
#' @section Display:
#' \describe{
#'  \item{\code{show(object)}:}{Display a peptideSet object.}
#'  \item{\code{summary(object)}:}{Summarize a peptideSet object.}
#' }
#'
#' @exportMethod "show"
#' @exportMethod "summary"
#'
#' @section Subset:
#' \describe{
#'  \item{\code{x[i, j]}:}{Subset x by peptides (i), or samples (j).}
#'  \item{\code{subset(x, subset, drop=FALSE)}:}{Subset x given an expression 'subset'.}
#' }
#'
#' @exportMethod "["
#' @exportMethod "subset"
#'
#' @importMethodsFrom IRanges lapply ranges ranges<- values values<- width cbind
#' rbind
NULL


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

# Subset by peptide/sample
setMethod("[", signature("peptideSet", i = "ANY", j = "ANY"),
          function(x, i, j, ..., drop = FALSE) {
            if (!missing(i)) {
              sdata <- exprs(x)[i, j, drop = drop]
              featureRange <- ranges(x)[i, ]
            } else {
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
		function(x, value)
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
				peptides <- values(ranges(x))[, type]
			}
      else {
				warning("'",type, "' is not a valid type! The accepted types are: ",
                paste(validTypes, collapse =", "),".")
			}
      return(peptides)
		})

setGeneric("peptide<-", function(object, value) standardGeneric("peptide<-"))
setReplaceMethod("peptide", signature("peptideSet", "character"), function(object, value){
  values(ranges(object))[, "peptide"] <- value
  return(object)
})

setGeneric("featureID", function(x, ...) standardGeneric("featureID"))
setMethod("featureID", "peptideSet",
		function(x, type=NULL){
			if (is.null(type))
			{
				values(ranges(x))[["featureID"]]
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
    exprs<-.bgCorrect.pSet(x)
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
          signature = signature(object = "GRanges"),
          definition = function(object){
            if(is.null(object$clade)){
              stop("The object does not have clade information!")
            }
            cladeList <- unique(unlist(
              strsplit(levels(as.factor(object$clade)), ","))) #List of all possible clades
            len <- length(object)
            retMatrix <- matrix(FALSE, nrow = len, ncol = length(cladeList))
            pepClades <- strsplit(as.character(object$clade), split = ",") #clades for each peptide
            for(pepIdx in 1:len){
              tmpList <- cladeList %in% pepClades[[pepIdx]]
              retMatrix[pepIdx, ] <- tmpList
            }
            rownames(retMatrix) <- names(object)
            colnames(retMatrix) <- cladeList
            return(retMatrix)
          })


setMethod("clade",
                signature = signature(object="peptideSet"),
                definition = function(object){ clade(object@featureRange)
                })

setGeneric("pepZscore", function(object) standardGeneric("pepZscore"))
setMethod("pepZscore", signature("GRanges"), function(object){
  vals <- as.data.frame(values(object))
  zs <- c("z1", "z2", "z3", "z4", "z5")
  zIns <- zs[zs %in% colnames(vals)]
  return(vals[,zIns])
})

setMethod("pepZscore", signature("peptideSet"), function(object){
  pepZscore(ranges(object))
})

setGeneric("pepZscore<-", function(object, value) standardGeneric("pepZscore<-"))
setReplaceMethod("pepZscore", signature("GRanges", "data.frame"), function(object, value){
  zs <- c("z1", "z2", "z3", "z4", "z5")
  if(!all(zs %in% colnames(value))){
    stop("The given data.frame does not contain the required colum names: 'z1', 'z2', 'z3', 'z4', 'z5'")
  }
  for(z in zs){
    object[[z]] <- value[[z]]
  }
  return(object)
})


setReplaceMethod("pepZscore", signature("peptideSet", "data.frame"), function(object, value){
  pepZscore(ranges(object)) <- value
  return(object)
})

setMethod("colnames", "peptideSet", function(x){ colnames(ranges(x)) })

# setReplaceMethod("rownames", signature("peptideSet", "character"),
#                  function(x, value){
#                    rownames(ranges(x)) <- value
#                    rownames(exprs(x)) <- value
#                    return(x)
#                    })
