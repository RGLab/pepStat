# This is a list of functions are simply used internaly, for now.
.bgCorrect.pSet <- function(pSet,verbose=FALSE){
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

#  .reduce2hotspots<-function(pSet,ranges,summary="median"){
#      FUN <- match.fun(summary)
#      y<-.bgCorrect.pSet(pSet)
#      sy<-sapply(1:nrow(hotspots),function(i,pSet,y,hotspots){apply(y[position(pSet)>start(hotspots[i,]) & position(pSet)<end(hotspots[i,]),],2,FUN)},pSet,y,hotspots)
#      # do.call(cbind,sy)
#      colnames(sy)<-rownames(ranges)
#      sy
#  }
#

# Checks that a peptideSet is valid
#   Order of samples
#   Order of peptides
.check_peptideSet <- function(peptideSet){
  if (class(peptideSet)!="peptideSet"){
    stop("peptideSet must be an object of class peptideSet")
  }
  if(any(colnames(exprs(peptideSet)) != rownames(pData(peptideSet)))){
    stop("The samples in the phenoData and assayData slots are different")
  }
  if(any(rownames(exprs(peptideSet)) != rownames(ranges(peptideSet)))){
    stop("The features in the featureRange and assayData slots are different")
  }
}

# Replace rownames in both featureRange and exprs slot
#   setReplaceMethod("rownames", ....) fails
.set_rownames <- function(x, value){
  names(ranges(x)) <- value
  rownames(exprs(x)) <- value
  return(x)
}


.makeZpepMatrix <- function(Sequence){
  Sequence = toupper(Sequence)
  let = unique(unlist(strsplit(Sequence, "")))
  AA = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
         "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  if(any(!(let %in% AA)))
    stop("Invalid sequences in Z-scale matrix construction")

  z = c(0.24, 0.84, 3.98, 3.11, -4.22,
        2.05, 2.47, -3.89, 2.29, -4.28,
        -2.85, 3.05, -1.66, 1.75, 3.52,
        2.39, 0.75, -2.59, -4.36, -2.54,
        -2.32, -1.67, 0.93, 0.26, 1.94,
        -4.06, 1.95, -1.73, 0.89, -1.3,
        -0.22, 1.62, 0.27, 0.5, 2.5,
        -1.07, -2.18, -2.64, 3.94, 2.44,
        0.6, 3.71, 1.93, -0.11, 1.06,
        0.36, 0.26, -1.71, -2.49, -1.49,
        0.47, 1.04, 1.84, -1.44, -3.5,
        1.15, -1.12, -1.54, 0.59, 0.43,
        -0.14, 0.18, -2.46, -3.04, 0.54,
        -0.82, 3.9, -0.84, 1.49, -0.72,
        1.94, -1.15, 0.7, -1.34, 1.99,
        -1.39, -1.46, -0.85, 3.44, 0.04,
        1.3, -2.65, 0.75, -0.25, -0.62,
        -0.38, 0.09, 0.26, 0.31, 0.84,
        -0.98, 1.61, 2, 0.66, -.17,
        0.67, -0.4, -0.02, -1.59, -1.47)
  dim(z) = c(20, 5)
  colnames(z) = paste0("z", 1:5)

  rownames(z) = AA

  Z = t(sapply(Sequence, .computeZpep, ztable = z))
  colnames(Z) = paste0("z", 1:5)
  rownames(Z) = Sequence
  Z
}

.computeZpep <- function(AAstring, ztable){
  if(AAstring == c("empty"))
    return(rep(0, 5))
  t = unlist(strsplit(AAstring, split = ""))
  colSums(ztable[t,])
}
