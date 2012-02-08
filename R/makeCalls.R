makeCalls<-function(peptideSet, cutoff=.1, method="local", freq=TRUE, group=NULL, verbose=FALSE)
{
  if(class(peptideSet)!="peptideSet")
  {
    stop("peptideSet must be an object of class peptideSet")
  }
  if(preproc(peptideSet@experimentData)$transformation!="log")
  {
    warning("The probe measurements should be log transformed!")
  }
  if(preproc(peptideSet@experimentData)$normalization=="none")
  {
    warning("You should probably normalize your data before using this function")
  }
  
  I<-.bgCorrect.pSet(peptideSet,verbose=verbose)
  
  if(method=="local")
  {
    Calls<-apply(I,2,.findFDR,cutoff,position)
  }
  else if(method=="global")
  {
    seqY<-seq(range(abs(I))[1],range(abs(I))[2],.1)
    FDR<-sapply(seqY,function(x,I){sum(I< -x)/sum(I>x)},as.double(I))
    Dmin<-seqY[which.min(abs(FDR-cutoff))]
    # print(cbind(FDR,seqY))
    Calls<-I>Dmin
  }
  else if(method=="absolute")
  {
    Calls<-I>cutoff
  }
  	  
#  browser()
  
  if(!is.null(group)&& freq) 
  {
	  	#parse the grouping variable 
		groupBy<-.parseCond(group)
		pd<-pData(peptideSet)
		#generate the factor list based on multipe grouping vairable
		factors<-lapply(groupBy,function(x){
					eval(substitute(pd$v,list(v=x)))
				})
		#split the ptid into groups
		ptidGroups<-split(pd$ptid,factors)
		#apply the rowMeans to each group
		res<-lapply(ptidGroups,function(curPtid){
#					browser()
					curCall<-Calls[,curPtid]
					rowMeans(curCall)
				})
	
	  return(res)
  }
  else if(freq)
  {
	  return(rowMeans(Calls))
  }
  Calls
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

.findFDR<-function(I,cutoff,position)
{
#	browser()
  seqY<-seq(range(abs(I))[1],range(abs(I))[2],.1)
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