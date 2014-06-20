##
#' Create a peptide collection
#'
#' Constructor to create peptide collection to be used in
#' \code{summarizePeptides}.
#'
#' @param position A \code{data.frame} or \code{GRanges} object. If a
#' \code{data.frame} is provided, it should contain 'start' and 'end' or 'width'
#' columns as well as a peptide column. If position is a \code{GRanges} object,
#' then it must either have peptide as names or contain a peptide metadata
#' column.
#'
#' @details
#'   \code{position} can have additional columns. These columns will be kept in
#'   the resulting peptide collection. This is especially useful to include
#'   clades and grouping parameters for the \code{makeCalls} function.
#'
#'   If the input contains all the z-scores (z1 to z5), then they will not be
#'   re-calculated. If some (but not all) z-scores are missing, a warning
#'   message will be sent and the z-scores are re-calculated.
#'
#'
#' @author Renan Sauteraud
#' @seealso \code{\link{GRanges}}
#'
#' @examples
#' #construct data.frame object
#'    AA <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P",
#'    "Q","R", "S", "T", "V", "W", "Y")
#'    starts <- seq(1, 30, 3)
#'    ends <- starts + 14
#'    peptides <- sapply(1:10, function(x) {
#'      paste0(AA[floor(runif(15, 1, 20))], collapse = "")
#'    })
#'    data <- data.frame(start = starts, end = ends, peptide = peptides)
#' #from data.frame
#'    new_pep <- create_db(data)
#' #from GRanges
#'    new_pep <- create_db(new_pep)
#'
#'
#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importMethodsFrom GenomicRanges mcols mcols<- names
#' @importClassesFrom IRanges IRanges
#' @importClassesFrom GenomicRanges GRanges
#'
create_db <- function(position){
  if(!class(position) %in% c("data.frame", "GRanges")){
    stop(paste("`position' must be of class 'GRanges' or data.frame'. Given:",
               class(position)))
  }else if(class(position) == "data.frame"){
    if(!"peptide" %in% colnames(position)){
      stop("The given data.frame must contain a 'peptide' column.")
    }
    cn <- colnames(position)
    extracols <- cn[!cn %in% c("start", "end", "width")]
    ir <- IRanges(start = position$start, end = position$end, width = position$width,
                  names = position$peptide)
    gr <- GRanges(ranges = ir, seqnames = " ")
    names(gr) <- position$peptide
    mcols(gr) <- position[, extracols, drop = FALSE]
  } else{
    if(is.null(names(position))){
      if(!"peptide" %in% colnames(mcols(position))){
        stop("The given GRanges object must contain a 'peptide' column.")
      } else{
        names(position) <- position$peptide
      }
    } else{
      if(class(names(position)) != "character"){
        stop("The names of the GRanges object should be the peptide sequences")
      }
    }
    if(any(names(position) != position$peptide)){
      stop("The names of the GRanges object should match the 'peptide' column.")
    }
    gr <- position
  }
  # Calculating z-scores
  if(!all(paste0("z", 1:5) %in% colnames(mcols(gr)))){
    existing_zs <- colnames(mcols(gr))[colnames(mcols(gr))%in% paste0("z", 1:5)]
    if(length(existing_zs) > 0){
      warning("Removing existing z-scores and recalculating new ones.")
      mcols(gr) <- mcols(gr)[!colnames(mcols(gr))%in% existing_zs]
    }
    mcols(gr) <- data.frame(mcols(gr), .makeZpepMatrix(gr$peptide))
  }
  return(gr)
}
# Test
# gr <- create_db(data)
# gr1 <- create_db(gr)
# gr2 <- gr; mcols(gr2) <- mcols(gr2)[, !colnames(mcols(gr2)) %in% paste0("z", 1:5), drop = FALSE]
# gr2 <- create_db(gr2)
# gr3 <- gr; gr3$z1 <- NULL
# gr3 <- create_db(gr3)
# d.f, complete GR, GR w/o zs, GR w/ missing z
