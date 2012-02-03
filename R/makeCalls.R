makeCalls<-function(peptideSet, cutoff=.1, method="local", average=TRUE, group=NULL,verbose=FALSE)
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
  
  pData(peptideSet)$ptid
  
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