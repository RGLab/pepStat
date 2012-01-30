makeCalls<-function(peptideSet, pos=NULL, cName=NULL, cutoff=.1, verbose=FALSE, method="local",paired=FALSE)
{
  if(class(peptideSet)!="peptideSet")
  {
    stop("peptideSet must be an object of class peptideSet")
  }
  if(preproc(peptideSet@experimentData)$transformation!="log" & preproc(peptideSet@experimentData)$transformation!="vsn")
  {
    stop("The probe measurements need to be log/vsn transformed!")
  }
  if(preproc(peptideSet@experimentData)$normalization=="none")
  {
    warning("You should probably normalize your data before using this function")
  }
  
  #whether the Control array is available
  y<-exprs(peptideSet)
  
  if(!is.null(cName))
  {
    sNames<-sampleNames(peptideSet)
    if(length(grep(cName,sNames))==0)
    {
      stop("The variable 'cName' must be a subset of the control sample name")
    }
    if(verbose)
    {
      cat("You are using: ",sNames[grep(cName,sNames)], " as the control\n" )
    }
    if(!paired)
    {
      C<-rowMeans(as.matrix(y[,grep(cName,sNames)]))
    }
    else
    {
      C<-as.matrix(y[,grep(cName,sNames)])
    }
    I<-as.matrix(y[,-grep(cName,sNames)])-C
  }
  else
  {
    I<-as.matrix(y)
  }

#browser()
  if(is.null(pos))
  {
    nProbes<-nrow(I)
    position<-1:nProbes
  }
  else
  {
	  #TODO:remove pos from argument list ?
#    pos<-pos[!duplicated(pos),]
#    match.seq<-match(peptide(peptideSet),pos$Sequence)
#    position<-pos$Position[match.seq]
#    df<-as.data.frame(I)
#  	tmp<-split(df,as.factor(position))
#    I<-t(sapply(tmp,function(x){apply(x,2,"mean")}))
#    position<-sapply(split(position,as.factor(position)),unique)
  }

  
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