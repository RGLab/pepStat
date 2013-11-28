#' peptideSet class
#' 
#' This class gathers all information from gpr files, annotation data and sequence data
#' 
#' @section Slots:
#' \describe{
#'  \item{featureRange:}{ A \code{RangedData object}. The ranges and sequences of 
#' the peptides and their associated annotation.}
#'  \item{phenoData:}{ An \code{AnnotatedDataFrame}. Annotation for the samples.}
#' }
#'
#' @seealso \code{\link{ExpressionSet}}, \code{\link{peptideSet-methods}}
#'
#' @importClassesFrom Biobase ExpressionSet
#' @importClassesFrom IRanges RangedData
#' @name peptideSet
#' @rdname peptideSet
#' @aliases peptideSet-class
#' @exportClass peptideSet
#' @author Greg Imholte
#' 
setClass("peptideSet",
    contains=c("ExpressionSet"),
    representation(featureRange="RangedData")
)

