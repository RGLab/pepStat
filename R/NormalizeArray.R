NormalizeArray<-function(peptideSet, robust=TRUE, standard=FALSE, method="ZpepQuad",
    centered=TRUE, verbose=FALSE)
{ 

  ### Sanity checks 
  if(class(peptideSet)!="peptideSet")
  {
    stop("peptideSet must be an object of class peptideSet (e.g. returned by makePeptideSet)!")
  }
  if(preproc(peptideSet@experimentData)$transformation!="log")
  {
    warning("The peptideSet measurements may not be log transformed!")
  }

  data<-exprs(peptideSet)  
  
  sNames<-sampleNames(peptideSet)
  # Initialize Zpep
  Zpep<-matrix(0,0,0)

  if(method=="constant")
  {
    methodNum<-1
  }
  else if(method == "binned")
  {
    methodNum=2
    numZpep = 0
  }
  else if(method == "Zpep")
  {
  	methodNum = 3
	X<-as.data.frame(ranges(peptideSet)@values[,grepl("z[1-9]",colnames(ranges(peptideSet)))][[1]])
	Zpep<-as.matrix(X)
  }
  else if(method == "ZpepQuad")
  {
	  methodNum = 3
	  X<-as.data.frame(ranges(peptideSet)@values[,grepl("z[1-9]",colnames(ranges(peptideSet)))][[1]])
	  Zpep<-as.matrix(cbind(X, X^2))
  }
  numZpep = 0
  if(!is.null(Zpep))
  {
  	if(!is.null(dim(Zpep)))
  	{
  		numZpep = ncol(Zpep)
  	}
  	else
  	{
  		numZpep = 1
  	}
  }
  
  y<-data
  nProbes<-nrow(y)
  nArrays<-ncol(y)

  seq<-peptide(peptideSet)

  # Extract the alphabet from the sequences
  alphabet<-sort(unique(AllChar<-unlist(sapply(seq, function(x){strsplit(x,"")}))))
  
  # Calling C code
  obj<-.C("NormalizeProbes",
  as.character(seq),
  as.double(y),
  yNormalized=as.double(rep(0,nArrays*nProbes)),
  as.integer(nProbes),
  as.integer(nArrays),
  as.integer(methodNum),
  as.integer(robust),
  adjRSquare=double(nArrays),
  RSquare=double(nArrays),
  BIC=double(nArrays),
  center=double(nArrays),
  beta=double(100*nArrays),
  betaLength=integer(1),
  as.character(alphabet),
  as.integer(length(alphabet)),
  as.integer(standard),
  as.integer(verbose),
  as.integer(numZpep),
  as.double(Zpep),
  package="pepStat")

  if(verbose)
  {
    cat("** Finished Normalizing ",nProbes," probes on ",nArrays," arrays **\n")
    cat("** Sample: R^2, BIC **\n")
    for(i in 1:length(obj$center))
    {
      cat("**", sNames[i],":", round(obj$RSquare[i],2),",", as.integer(obj$BIC[i])," **\n")
    }
  }
  # beta<-matrix(obj$beta[obj$beta!=0],length(alphabet)+1,ncol(y),byrow=TRUE)
  # save("beta",file="beta.rda")
#	  browser()
  normData<-matrix(obj$yNormalized, nProbes, length(obj$center))
  if(length(obj$center)==1)
  {
    # Reconvert to a matrix, in case there is only one array
    normData<-matrix(normData,nProbes,1)
  }

  
  if(!centered)
  {
    # Recenter the data
    normData<-as.matrix(normData)+rep(1,nProbes)%*%t(obj$center)
  }

  ### Set the normalized data as exprs
  exprs(peptideSet)<-normData
  colnames(exprs(peptideSet))<-sampleNames(peptideSet)
  ### Setting the normalization parameters
  normalization<-"Z-scale"
  if(robust)
  {
    normalization<-paste(normalization, "robust",sep=" ")
  }
  if(standard)
  {
    normalization<-paste(normalization, "standardized",sep=" ")
  }

  preproc(peptideSet@experimentData)$normalization<-normalization
  colnames(exprs(peptideSet))<-sampleNames(peptideSet)
  peptideSet
}
