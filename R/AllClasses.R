#' peptideSet class
#'
#' This class gathers all information from gpr files, annotation data and sequence data
#'
#' @section Slots:
#' \describe{
#'  \item{featureRange}{A \code{GRanges}object. The ranges and sequences of
#' the peptides and their associated annotation.}
#'  \item{phenoData}{An \code{AnnotatedDataFrame}. Annotation for the samples.}
#'  \item{assayData}{}
#'  \item{featureData}{}
#'  \item{annotation}{}
#'  \item{protocolData}{Slots inherited from \code{ExpressionSet}.}
#' }
#'
#' @details
#' See \code{?`peptideSet-methods`} for a list of accessors and method associated
#' with the class.
#'
#' @seealso \code{\link{ExpressionSet}}, \code{\link{peptideSet-methods}}
#'
#' @importFrom Biobase ExpressionSet
#' @importClassesFrom Biobase Versioned VersionedBiobase eSet ExpressionSet
#' @importClassesFrom GenomicRanges GRanges
#' @name peptideSet
#' @rdname peptideSet
#' @aliases peptideSet-class
#' @exportClass peptideSet
#' @author Greg Imholte
#'
setClass("peptideSet",
         contains=c("ExpressionSet"),
         representation(featureRange="GRanges")
)