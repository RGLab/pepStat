makeCalls <- function(peptideSet, cutoff=.1, method="absolute", freq=TRUE, group=NULL, verbose=FALSE)
{
	if (class(peptideSet)!="peptideSet") {
		stop("peptideSet must be an object of class peptideSet")
	}
	
	if (preproc(peptideSet@experimentData)$transformation!="log") {
		warning("The probe measurements should be log transformed.")
	}
	
	if (preproc(peptideSet@experimentData)$normalization=="none") {
		warning("You should probably normalize your data before using this function.")
	}
  
	I <- .bgCorrect.pSet(peptideSet, verbose=verbose)
  
	if (method == "FDR") {
		Calls<-.findFDR(I, cutoff, position)
		
		} else if(method == "absolute") {
		Calls <- I > cutoff
	}
  	  
  
	if (!is.null(group) && freq) {
		#parse the grouping variable 

		# Only select the Post and remove empty levels        
		t1 <- grepl("post", pData(peptideSet)$visit)
		pd <- drop.levels(pData(peptideSet)[t1, ])
        if(!group%in%colnames(pd)) {
            error("The grouping variable is not part of the pData object.")
        } else {
            factor<-as.factor(pd[,group])
        }
            
            
		if (nlevels(factor) > 1) {
			#split the ptid into groups
			ptidGroups <- split(pd$ptid,factor)
			#apply the rowMeans to each group
			res <- lapply(ptidGroups, function(curPtid,Calls,ptid){rowMeans(Calls[,ptid%in%curPtid])}, Calls, pd$ptid)
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

.findFDR <- function(I, cutoff, position)
{
    seqY <- seq(min(abs(I)), max(abs(I)),.05)
    # Split the data by unique positions
    tmp <- split(as.data.frame(I), position)
    # Compute the mean over unique positions
    D <- sapply(tmp,apply, 2, mean)
    
    FDR <- sapply(seqY, function(x, D){median(apply(D, 1, function(D,x){min(max(sum(D < -x)/sum(D > x), 0), 1)}, x), na.rm=TRUE)}, D)
    
    # Did not find anything below the cutoff or everything is NA
    if (all(round(FDR,2)>cutoff, na.rm=TRUE) | all(is.na(FDR))) {
        return(I > max(I)) # Return all FALSE
    } else {
        Dmin <- seqY[which.min(abs(FDR-cutoff))]
        return(I > Dmin)
    }
}