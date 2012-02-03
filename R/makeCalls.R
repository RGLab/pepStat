makeCalls<-function(peptideSet, cutoff=.1, method="local", freq=TRUE, group=TRUE, verbose=FALSE)
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
  
  y<-exprs(peptideSet)
  
  # if(!is.null(cName))
  # {
  #   sNames<-sampleNames(peptideSet)
  #   if(length(grep(cName,sNames))==0)
  #   {
  #     stop("The variable 'cName' must be a subset of the control sample name")
  #   }
  #   if(verbose)
  #   {
  #     cat("You are using: ",sNames[grep(cName,sNames)], " as the control\n" )
  #   }
  #   if(!paired)
  #   {
  #     C<-rowMeans(as.matrix(y[,grep(cName,sNames)]))
  #   }
  #   else
  #   {
  #     C<-as.matrix(y[,grep(cName,sNames)])
  #   }
  #   I<-as.matrix(y[,-grep(cName,sNames)])-C
  # }
  # else
  # {
  #   I<-as.matrix(y)
  # }
  
  ptid<-pData(peptideSet)$ptid
  t0<-grepl("[Pp][Rr][Ee]",pData(peptideSet)$visit)
  t1<-grepl("[Pp][Oo][Ss][Tt]",pData(peptideSet)$visit)
  ### Paired
  if(all.equal(sort(ptid[t0]),sort(ptid[t1])))
  {
	  if(verbose)
	  {
	  	cat("You have paired PRE/POST samples\n")
	  }
	  I<-as.matrix(y[,t1])-as.matrix(y[,t0])
  }
  else
  {
	  if(verbose)
	  {
	  	cat("You don't have paired PRE/POST samples\n")
	  }	  
	  I<-as.matrix(y[,t1])-as.matrix(rowMeans(y[,t0]))  	
  }
  
  colnames(I)<-ptid[t1]
  rownames(I)<-peptide(peptideSet)
  
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
  
  if(group && nlevels(group)>1 && freq)
  {
	  group<-as.factor(pData(peptideSet)$treatment)[t1]
	  freq<-as.data.frame(sapply(levels(group),function(x,Calls,group){rowMeans(Calls[,group==x])},Calls,group))
	  names(freq)<-levels(group)
	  return(freq)
  }
  else if(freq)
  {
	  return(rowMeans(Calls))
  }
  Calls
}

.findFDR<-function(I,cutoff,position)
{
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