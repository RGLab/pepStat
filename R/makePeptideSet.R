makePeptideSet<-function(files=NULL, path=NULL, use.flags=TRUE, rm.control.list=c("empty","none","JPT-","Ig","Cy","landmark"), norm.empty=TRUE, empty.control.list=c("empty","blank control"),bgCorrect.method="normexp", log=TRUE, verbose=FALSE)
{
	# There is some ambiguity with respect to what is Name and ID
	# ID -> peptide
	# Name -> annotation
	f <- function(x) as.numeric(x$Flags > -99)
	# Only read the red channel
	RG<-read.maimages(files=files,source="genepix", path=path, ext="gpr", columns=list(R="F635 Median",Rb="B635 Median"), wt.fun=f,verbose=verbose)

	RG<-backgroundCorrect(RG, method=bgCorrect.method, offset=1, verbose=verbose)

	myDesc <- new("MIAME")

	if(use.flags)
	{
		RG$E[RG$weights==0]<-NA
	}
	if(log)
	{
		RG$E<-log2(RG$E)
		preproc(myDesc)<-list(transformation="log", normalization="none")
	}
	else
	{
		preproc(myDesc)<-list(transformation="none", normalization="none")
	}

	if(norm.empty)
	{
		#Check both name and ID for the control list
		index<-RG$genes$Name%in%empty.control.list | RG$genes$ID%in%empty.control.list
		if(verbose)
		{
			cat("** Background corrected using ", sum(index), " empty splots **")
		}		
		mean.empty<-matrix(colMeans(as.matrix(RG$E[index,])),nrow=nrow(as.matrix(RG$E)),ncol=ncol(as.matrix(RG$E)),byrow=TRUE)
	}
	else
	{
		mean.empty<-rep(0,ncol(as.matrix(RG$E)))
	}

	if(!is.null(rm.control.list))
	{
		ind.keep<-lapply(rm.control.list,function(x,Name,ID){!grepl(x,Name) & !grepl(x,ID)},as.character(RG$genes$Name),as.character(RG$genes$ID))
		ind.keep<-do.call(cbind,ind.keep)
		ind.keep<-apply(ind.keep,1,all)
	}

	## See note above about Name and ID
	featureSequence<-as.character(RG$genes$ID)[ind.keep]
	featureID<-as.character(RG$genes$Name)[ind.keep]
	nPep<-length(which(ind.keep))
	#browser()
	new('peptideSet'
	,featureRange=RangedData(IRanges(rep(0,nPep),rep(0,nPep))
	,featureID
	,peptide=featureSequence
	)
	,exprs=as.matrix(RG$E-mean.empty)[ind.keep,]
	,experimentData=myDesc
	)
			
}

