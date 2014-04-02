inter_group <- function(peptideSet, calls){
  if(!(all(rownames(calls) == rownames(peptideSet)))){
    stop("The call result and peptideSet have different peptides")
  }
  data <- data.table(calls)
  data <- data[, position := position(peptideSet)]
  data <- data[, lapply(.SD, mean), by = "position"]
  data <- data[order(position)]
  return(data)
}

##
#' Plot frequency of response for each group
#' 
#' Plot an axis and the frequency of response of each group, averaged by peptides
#' at each position.
#' 
#' @param peptideSet A \code{peptideSet} object.
#' @param calls A \code{matrix}, as returned by the \code{makeCalls} function.
#' @param from A \code{numeric}, the start coordinate of the plot.
#' @param to A \code{numeric}, the end coordinate of the plot.
#' 
#' @details 
#' The peptideSet should be the one used in the function call to \code{makeCalls}
#' that generated the calls used. They should have identical peptides in the same
#' order.
#' The function requires the calls matrix to have been created using \code{makeCalls}
#' with freq set to TRUE and a grouping variable.
#' 
#' @seealso plot_clade
#' @author Renan Sauteraud
#' @export
#' 
plot_inter <- function(peptideSet, calls, from=0, to){
  data <- inter_group(peptideSet, calls)
  if(missing(to)){
    to <- max(data$position)
  } else if(to > max(data$position)){
    warning(paste("'to' is bigger than the sequence length. Set to", max(data$position)))
    to <- max(data$position)
  }
  if(from > to){
    stop(paste0("'from' (", from, ") is bigger than 'to' (", to, ")"))
  }
  pat <- ProteinAxisTrack(addNC = TRUE)
  dt <- DTrack(start=data$position, end=data$position, data=data[, 2:ncol(data), with=FALSE],
              groups =  colnames(data)[2:ncol(data)], name="Freq", legend=TRUE, type="l", background.title="darkgray")
  plotTracks(c(pat, dt), from = from, to = to)
}

##
#' Plot frequency of response for a single clade
#' 
#' Plot an axis and the frequency of response of a single selected clade
#' 
#' @param peptideSet A \code{peptideSet} object.
#' @param calls A \code{matrix}, as returned by the \code{makeCalls} function.
#' @param clade A \code{character}, the name of the clade to plot. It must be a 
#' valid column name of \code{clade(peptideSet)}.
#' @param from A \code{numeric}, the start coordinate of the plot.
#' @param to A \code{numeric}, the end coordinate of the plot.
#' 
#' @details 
#' The peptideSet should be the one used in the function call to \code{makeCalls}
#' that generated the calls used. They should have identical peptides in the same
#' order.
#' The function requires the calls matrix to have been created using \code{makeCalls}
#' with freq set to TRUE and a grouping variable.
#' 
#' @seealso plot_inter
#' @author Renan Sauteraud
#' @export
#' 
plot_clade <- function(peptideSet, calls, clade, from = 0, to=max(end(peptideSet))){
  CLADE <- clade
  data <- data.frame(calls)
  data$position <- pepStat::position(peptideSet)
  data$clade <- ranges(peptideSet)$clade
  data <- data[grep(CLADE, data$clade),]
  
#   if(missing(to)){
#     to <- max(data$position)
#   } else if(to > max(data$position)){
#     warning(paste("'to' is bigger than the sequence length. Set to", max(data$position)))
#     to <- max(data$position)
#   }
  if(from > to){
    stop(paste0("'from' (", from, ") is bigger than 'to' (", to, ")"))
  }
  pos <- data$position
  data <- data[, grep("^(position|clade)$", colnames(data), invert=TRUE)]
  pat <- ProteinAxisTrack(addNC = TRUE)
  print("HELLO")
 pdt <- DTrack(start=pos, end=pos, data=data)#, groups = colnames(data), 
                 #name=CLADE, legend=TRUE, type="l", background.title = "darkgray")
print(paste(clade, from, to))
#     plotTracks(c(pat, dt), from=from, to=to)
}