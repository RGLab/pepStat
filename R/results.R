#' Result table
#'
#' Tabulate the results of a peptide microarray analysis.
#'
#' @param peptideSet A \code{peptideSet} object.
#' @param calls A \code{matrix}, as returned by the \code{makeCalls} function.
#'
#' @details
#' The peptideSet should be the one used in the function call to \code{makeCalls}
#' that generated the calls used. They should have identical peptides.
#'
#' @return A \code{data.frame} with the peptides and some information from the
#' \code{peptideSet} as well as the frequency of binding for each group of the
#' calls.
#'
#' @export
#' @importMethodsFrom GenomicRanges as.data.frame
#'
#' @example examples/pipeline.R
restab <- function(peptideSet, calls){
  pep <- as.data.frame(ranges(peptideSet))
  pep$names <- rownames(pep)
  pep$position <- pepStat::position(peptideSet)
  cn <- c("names", "peptide", "position", "start", "end", "width", "clade")
  pep <- pep[, cn]
  calls <- data.frame(calls)
  calls$names <- rownames(calls)
  restab <- merge(pep, calls, by = "names")
  restab <- restab[order(restab$position),]
  rownames(restab) <- restab$names
  restab$names <- NULL
  return(restab)
}
