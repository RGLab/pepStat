### These are some functions that are mostly for internal use. We might expose them in the future
.reshape.pSet<-function(pSet)
{
	ind<-grepl("[Pp][Rr][Ee]",pData(pSet)$visit)
	y<-as.data.frame(exprs(pSet))[,ind,drop=FALSE]	
	names(y)<-pData(pSet)$ptid[ind]
	tmp1<-melt(y,value.name="Intensity",variable.name="ptid")
	tmp1<-cbind(tmp1,visit="pre")
	y<-as.data.frame(exprs(pSet))[,!ind,drop=FALSE]	
	names(y)<-pData(pSet)$ptid[!ind]	
	tmp2<-melt(y,value.name="Intensity",variable.name="ptid")
	tmp2<-cbind(tmp2,visit="post")
	y.long<-rbind(tmp1,tmp2)
}

.bgCorrect.pSet<-function(pSet,verbose=FALSE)
{
    y<-exprs(pSet)
    ptid<-pData(pSet)$ptid
    t0<-grepl("[Pp][Rr][Ee]",pData(pSet)$visit)
    t1<-grepl("[Pp][Oo][Ss][Tt]",pData(pSet)$visit)
    ### Paired
#	browser()
	if(length(ptid[t0])==0||length(ptid[t1])==0)
	{
		I<-as.matrix(y[,t1])
		
	}else
	{
		if(isTRUE(all.equal(sort(ptid[t0]),sort(ptid[t1]))))
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
			I<-as.matrix(y[,t1])-rowMeans(y[,t0,drop=FALSE],na.rm=TRUE)#the vector to be subtracted from matrix need to be the same length as nrow of the matrix  	
		}
	}
	colnames(I)<-ptid[t1]
    rownames(I)<-peptide(pSet)
	I	
}

.reduce2hotspots<-function(pSet,ranges,summary="median")
{
    FUN <- match.fun(summary)
	y<-pepStat:::.bgCorrect.pSet(pSet)
	sy<-sapply(1:nrow(hotspots),function(i,pSet,y,hotspots){apply(y[position(pSet)>start(hotspots[i,]) & position(pSet)<end(hotspots[i,]),],2,FUN)},pSet,y,hotspots)
	# do.call(cbind,sy)
	colnames(sy)<-rownames(ranges)
	sy
}

.impute.pSet<-function(pSet,cv=FALSE,k=20,verbose=T)
{
### ute missing data based on 	
    y<-exprs(pSet)
    ptid<-pData(pSet)$ptid
    t0<-grepl("[Pp][Rr][Ee]",pData(pSet)$visit)
    t1<-grepl("[Pp][Oo][Ss][Tt]",pData(pSet)$visit)
	if(length(ptid[t0])==0||length(ptid[t1])==0)
	{
		I<-as.matrix(y[,t1])
		
	}else
	{
		if(isTRUE(all.equal(sort(ptid[t0]),sort(ptid[t1]))))
		{
			C<-as.matrix(y[,t0,drop=FALSE])
		}
		else
		{
			C<-rowMeans(y[,t0],na.rm=TRUE)#the vector to be subtracted from matrix need to be the same length as nrow of the matrix  	
		}
	}
	I<-as.matrix(y[,t1])-C	
	
	if(cv)
	{
		k<-cv.kNNute(I,k.max=k)$k
		if(verbose)
		{
			cat("k.best=",k,"\n")
		}
	}
	I.new<-kNNImpute(I,k=k,verbose=verbose)$x

	### Only impute the post
	y[,t1]<-I.new+C
	
	exprs(pSet)<-y
	pSet
}
