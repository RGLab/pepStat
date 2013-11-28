##
#' Create a peptide collection
#'
#' Constructor to create peptide collection such as the datasets available in PEP.db.
#'
#' @param rd A \code{RangedData} object. The object should have a peptide column.
#'
#' @details
#'   rd can have additional columns. These columns will be kept in the peptide collection.
#' 
#' @seealso \code{\link{RangedData}}
#' 
#' 
#'
#' @export
#' @importClassesFrom IRanges RangedData
##

create_db <- function(rd){
  if(class(rd) != "RangedData"){
    stop("`rd' must be an object of class 'RangedData'. Given: ",class(rd))
  }
  if(is.null(rd$peptide)){
    stop("The RangedData should have a peptide column")
  }
  zs <- makeZpepMatrix(as.character(rd$peptide))
  rownames(zs) <- NULL
  pepZscore(rd) <- as.data.frame(zs)
  rownames(rd) <- as.character(rd$peptide)
  return(rd)
}
