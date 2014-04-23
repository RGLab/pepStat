##
#' Create a peptide collection
#'
#' Constructor to create peptide collection such as the datasets available in PEP.db.
#'
#' @param rd A \code{RangedData} object. The object should have a peptide column.
#'
#' @details
#'   rd can have additional columns. These columns will be kept in the peptide
#'   collection.
#'
#' @seealso \code{\link{RangedData}}
#'
#' @examples
#' #construct RangedData object
#'    library(IRanges)
#'    AA <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P",
#'    "Q","R", "S", "T", "V", "W", "Y")
#'    starts <- seq(1, 30, 3)
#'    ends <- starts + 14
#'    peptides <- sapply(1:10, function(x) {
#'      paste0(AA[floor(runif(15, 1, 20))], collapse = "")
#'    })
#'    data <- data.frame(start = starts, end = ends, peptide = peptides)
#'    newRD <- RangedData(data)
#' #create_db
#'    new_pep <- create_db(newRD)
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
  if("clade" %in% colnames(rd) && class(rd$clade) != "character"){
    rd$clade <- as.character(rd$clade)
  }
  return(rd)
}
