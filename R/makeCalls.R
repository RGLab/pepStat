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
		Calls<-apply(I, 2, .findFDR, cutoff, position)
		
		} else if(method == "absolute") {
		Calls <- I>cutoff
	}
  	  
  
	if (!is.null(group) && freq) {
		#parse the grouping variable 
		groupBy <- .parseCond(group)
		t1 <- grepl("post", pData(peptideSet)$visit)
		#Only select the Post and remove empty levels
		pd <- drop.levels(pData(peptideSet)[t1, ])
		
		#generate the factor list based on multipe grouping vairable
		factors <- lapply(groupBy,function(x,pd){eval(substitute(pd$v,list(v=x)))},pd)[[1]]
		
		if (nlevels(factors) > 1) {
			#split the ptid into groups
			ptidGroups <- split(pd$ptid,factors)
			#apply the rowMeans to each group
			res <- lapply(ptidGroups, function(curPtid,Calls,ptid){rowMeans(Calls[,ptid%in%curPtid])}, Calls, pd$ptid)
			res <- do.call(cbind, res)
		} else {
			return(rowMeans(Calls)*100)
		}
		return(res*100)
	} else if(freq) {
        return(rowMeans(Calls)*100)
	}
}

.parseCond <-
function(model)
{
    ## WAS: model <- eval(parse(text = paste("~", deparse(model))))[[2]]
    ## but that's not good (PR#7395)
    model <- substitute(~m, list(m = model))[[2]]
    model.vars <- list()
    while (length(model) == 3 && (model[[1]] == as.name("*")
    || model[[1]] == as.name("+"))) {
        model.vars <- c(model.vars, model[[3]])
        model <- model[[2]]
    }
    rev(c(model.vars, model))
}

.getFDR <- function(y, t)
{
    min(max(sum(y < -t)/sum(y > t),0),1)
}

.findFDR<-function(I,cutoff,position)
{
    #	browser()
    seqY<-seq(range(abs(I))[1],range(abs(I))[2],.05)
    tmp<-split(as.data.frame(I),position)
    D<-unlist(sapply(tmp,function(x){rep(apply(x,2,mean),nrow(x))}))
    FDR<-sapply(seqY,function(x,D){sum(D< -x)/sum(D>x)},D)
    # print(cbind(seqY,FDR))
    if(!any(round(FDR,2)<=cutoff))
    {
        return(rep(FALSE,length(I)))
    }
    else
    {
        Dmin<-seqY[which.min(abs(FDR-cutoff))]
        return(I>Dmin)
    }
}