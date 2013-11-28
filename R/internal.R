# This is a list of functions are simply used internaly, for now.
.bgCorrect.pSet <- function(pSet,verbose=FALSE)
{
    y <- exprs(pSet)
    ptid <- pData(pSet)$ptid
    t0 <- grepl("[Pp][Rr][Ee]", pData(pSet)$visit)
    t1 <- grepl("[Pp][Oo][Ss][Tt]", pData(pSet)$visit)

    ### Paired
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
    rownames(I) <- peptide(pSet)
	I
}
.bgCorrect2.pSet <- function(pSet, verbose=FALSE){
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
  rownames(I) <- peptide(pSet)
  I
}
  


 .reduce2hotspots<-function(pSet,ranges,summary="median")
 {
     FUN <- match.fun(summary)
     y<-.bgCorrect.pSet(pSet)
     sy<-sapply(1:nrow(hotspots),function(i,pSet,y,hotspots){apply(y[position(pSet)>start(hotspots[i,]) & position(pSet)<end(hotspots[i,]),],2,FUN)},pSet,y,hotspots)
     # do.call(cbind,sy)
     colnames(sy)<-rownames(ranges)
     sy
 }
 
