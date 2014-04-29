#' Make antibody binding positivity calls
#'
#' After normalization and data smoothing, this last step makes the call for each
#' peptide of the peptideSet after baseline correcting the peptide intenstities.
#'
#' @param peptideSet A \code{peptideSet} object. The peptides, after normalization
#' and possibly data smoothing.
#' @param cutoff A \code{numeric}. If FDR, the FDR threshold. Otherwise, a cutoff
#' for the background corrected intensities.
#' @param method A \code{character}. The method used to make positivity calls.
#' "absolute" and "FDR" are available. See details below.
#' @param freq A \code{logical}. If set to TRUE, return the percentage of slides
#' calling a peptide positive. Otherwise, return a \code{logical} indicating binding
#' events.
#' @param group A \code{character}. Only used when freq is set to TRUE. A character indicating
#' a variable by which to group slides. If non-null the percentage is calculated by group.
#' @param verbose A \code{logical}. If set to TRUE, progress information will be displayed.
#'
#' @details
#' This function requires specific variables ptid and visit in pData(peptideSet).
#' The variable \code{ptid} should indicate subjects, and the variable \code{visit}
#' should be a factor with levels pre and post.
#'
#' If slides are paired for subjects, intensities corresponding to post-visit are
#' substracted from pre. If slides are not paired, slides with pre have intensities
#' averaged by peptides, and averaged peptide intensities are subtracted from slides
#' that have entry post. Calls are made on these baseline corrected intensities.
#'
#' When method = FDR, a left-tail method is used to generate a threshold controlling
#' the False Discovery Rate at level \code{cutoff}. When method = absolute, Intensities
#' exceeding the threshold are labelled as positive.
#'
#' When freq = TRUE a group variable may be specified. The argument group indicates
#' the name of a variable in pData(peptideSet) by which positive calls should be grouped.
#' The call frequency for each peptide is calculated within groups.
#'
#' @return If freq = TRUE, a \code{numeric} \code{matrix} with peptides as rows and
#' groups as columns where the values are the frequency of response in the group. If
#' freq = FALSE, a \code{logical} \code{matrix} indicating binding events for each
#' peptide in each subject.
#'
#' @rdname makeCalls
#' @aliases makeCalls
#' @author Greg Imholte
#' @export
#' @example examples/pipeline.R
makeCalls <- function(peptideSet, cutoff=1.2, method="absolute", freq=TRUE, group=NULL, verbose=FALSE){
	if (class(peptideSet)!="peptideSet") {
		stop("peptideSet must be an object of class peptideSet")
	}

	if (preproc(peptideSet@experimentData)$transformation!="log") {
		warning("The probe measurements should be log transformed.")
	}

	if (preproc(peptideSet@experimentData)$normalization=="none") {
		warning("You should probably normalize your data before using this function.")
	}

	I <- baselineCorrect.pSet(peptideSet, verbose=verbose)

	if (method == "FDR") {
		Calls<-.findFDR(I, cutoff, position(peptideSet), verbose=verbose)
	} else if(method == "absolute") {
		Calls <- I > cutoff
	}

	if (!is.null(group) && freq) {
		#parse the grouping variable

		# Only select the Post and remove empty levels
		t1 <- grepl("post", tolower(pData(peptideSet)$visit))
		pd <- pData(peptideSet)[t1, ]

    if(!group%in%colnames(pd)){
      stop("The grouping variable is not part of the pData object.")
      } else {
        factor<-factor(pd[,group])
      }
		if(nlevels(factor) > 1){
			#split the ptid into groups
			ptidGroups <- split(pd$ptid,factor)
			#apply the rowMeans to each group
			res <- lapply(ptidGroups, function(curPtid,Calls,ptid){rowMeans(Calls[,ptid%in%curPtid, drop=FALSE])}, Calls, pd$ptid)
			res <- do.call(cbind, res)
		} else {
			return(rowMeans(Calls)*100)
		}
    return(res*100)
	} else if(freq) {
        return(rowMeans(Calls)*100)
    } else {
        return(Calls)
    }
}

.findFDR <- function(I, cutoff, position, verbose=FALSE){
    seqY <- seq(min(abs(I)), max(abs(I)),.05)
    # Split the data by unique positions
    tmp <- split(as.data.frame(I), position)
    # Compute the mean over unique positions
    D <- sapply(tmp,apply, 2, mean)

    #Calculate F(T)
    FDR <- sapply(seqY, function(x, D){median(
                  #Fs(T)
                  apply(D, 1, function(D,x){min(sum(D < -x)/sum(D > x), 1)}, x), na.rm=TRUE)}, D)

    # Did not find anything below the cutoff or everything is NA
    if (all(round(FDR,2)>cutoff, na.rm=TRUE) | all(is.na(FDR))) {
        return(I > max(I)) # Return all FALSE
    } else {
        Dmin <- seqY[which.min(abs(FDR-cutoff))]
        if(verbose){
          cat("The selected threshold T is ", Dmin, "\n")
        }
        return(I > Dmin)
    }
}


#' Substract baseline intensities
#'
#' Correct intensities by substracting PRE visit sample intensities.
#'
#' @param pSet A \code{peptideSet} with sample PRE and POST visits.
#' @param verbose A \code{logical}. If TRUE, information regarding the
#' pairedness of the data will be displayed.
#'
#' @return A \code{matrix} of the baseline corrected intensities, with as many
#' columns as there are samples POST visit
#'
#' @details
#' If samples are not PAIRED (One PRE and POST for each ptid), then the average
#' expression of all PRE visit samples is substracted from each sample.
#'
#' @author Raphael Gottardo, Gregory Imholte
#'
baselineCorrect.pSet <- function(pSet, verbose=FALSE){
  y <- exprs(pSet)
  ptid <- pData(pSet)$ptid
  t0 <- grep("pre", tolower(pData(pSet)$visit))
  t1 <- grep("post", tolower(pData(pSet)$visit))
  ### Paired?
  if (length(ptid[t0])==0||length(ptid[t1])==0) {
    I<-as.matrix(y[,t1])
  } else {
    if (isTRUE(all.equal(sort(ptid[t0]), sort(ptid[t1])))) {
      if (verbose) {
        message("You have paired PRE/POST samples\n")
      }
      I <- as.matrix(y[,t1])-as.matrix(y[,t0])
    } else {
      if(verbose) {
        message("You don't have paired PRE/POST samples\n")
      }
      I <- as.matrix(y[,t1])-rowMeans(y[, t0, drop=FALSE], na.rm=TRUE)#the vector to be subtracted from matrix need to be the same length as nrow of the matrix
    }
  }
  colnames(I) <- ptid[t1]
  if(isTRUE(preproc(pSet)$split.by.clade)){
    rownames(I) <- rownames(pSet)
  } else{
    rownames(I) <- peptide(pSet)
  }
  return(I)
}

#' Substract baseline intensities
#'
#' Correct intensities by substracting PRE visit sample intensities.
#'
#' @param pSet A \code{peptideSet} with sample PRE and POST visits.
#' @param verbose A \code{logical}. If TRUE, information regarding the
#' pairedness of the data will be displayed.
#'
#' @return A \code{matrix} of the baseline corrected intensities, with as many
#' columns as there are samples POST visit
#'
#' @details
#' The function will try to pair as many sample as possible. The remaining subjects
#' with a POST and no PRE will use the average expression of all baseline samples.
#' Subjects with baseline only will not be represented in the resulting matrix.
#'
#' @author Renan Sauteraud
baseline_correct <- function(pSet, verbose=FALSE){
  exprs <- exprs(pSet)
  pd <- pData(pSet)
  t0 <- grep("pre", tolower(pd$visit))
  t1 <- grep("post", tolower(pd$visit))
  ## All paired
  if (isTRUE(all.equal(sort(pd$ptid[t0]), sort(pd$ptid[t1])))){
    if (verbose) {
      message("You have paired PRE/POST samples\n")
    }
    I <- as.matrix(exprs[,t1])-as.matrix(exprs[,t0])
  } else {
    ## Deal with paired samples first
    paired <- pd$ptid[t1][pd$ptid[t1] %in% pd$ptid[t0]]
    pre_paired <- rownames(pd[ pd$ptid %in% paired & tolower(pd$visit)=="pre",])
    post_paired <- rownames(pd[ pd$ptid %in% paired & tolower(pd$visit)=="post",])
    I <- exprs[, post_paired] - exprs[, pre_paired]
    colnames(I) <- paired
    ## Cbind the rest
    if (verbose) {
      cat(length(non_paired), "subjects don't have a baseline sample.\n")
    }
    non_paired <- unique(pd$ptid[!( pd$ptid %in% paired) & tolower(pd$visit)=="post"])
    if(length(non_paired) > 0){
      post_only <- rownames(pd[ pd$ptid %in% non_paired, ])
      I_no_pairs <- exprs[, post_only] - rowMeans(exprs[, t0, drop=FALSE], na.rm=TRUE)
      colnames(I_no_pairs) <- non_paired
      I <- cbind(I, I_no_pairs)
    }
  }
  rownames(I) <- peptide(pSet)
  return(I)
}
