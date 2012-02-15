.reshape.pSet<-function(pSet)
{
	ind<-grepl("[Pp][Rr][Ee]",pData(pSet)$visit)
	y<-as.data.frame(exprs(pSet))[,ind]	
	names(y)<-pData(pSet)$ptid[ind]
	tmp1<-melt(y,value.name="Intensity",variable.name="ptid")
	tmp1<-cbind(tmp1,visit="pre")
	y<-as.data.frame(exprs(pSet))[,!ind]	
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
			I<-as.matrix(y[,t1])-rowMeans(y[,t0],na.rm=TRUE)#the vector to be subtracted from matrix need to be the same length as nrow of the matrix  	
		}
	}
	colnames(I)<-ptid[t1]
    rownames(I)<-peptide(pSet)
	I	
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
			C<-as.matrix(y[,t0])
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

.plot.aggregate<-function(pSet,freq,anno.list,minimalist=TRUE,colorTrack,colorAnno,sizeTrack,sizeAnno,hotspots)
{
    t0<-grepl("[Pp][Rr][Ee]",pData(pSet)$visit)
    t1<-grepl("[Pp][Oo][Ss][Tt]",pData(pSet)$visit)
	
	
	freq.track<-makeGenericArray(as.matrix(freq),probeStart=position(pSet),dp=DisplayPars(size=sizeTrack, color = paste(colorTrack[1:2],95,sep=""), type="line",lwd=4,cex=.2,axisCex=.7,isLegend=is.matrix(freq)))
	y<-pepStat:::.bgCorrect.pSet(pSet,verbose=FALSE)
	y.smooth<-makeSmoothing(position(pSet), rowMeans(y), dp = DisplayPars(color = "#80808090", type="line", lwd=4))
	y.track<-makeGenericArray(as.matrix(y), probeStart=position(pSet), dp = DisplayPars(size=sizeTrack, color = "#00000080", type="line",lwd=1),trackOverlay=list(y.smooth))
	y.post<-exprs(pSet[,t1])
	y.post.smooth<-makeSmoothing(position(pSet), rowMeans(y.post), dp = DisplayPars(color = "#80808090", type="line",lwd=4))
	y.post.track<-makeGenericArray(as.matrix(y.post), probeStart=position(pSet), dp = DisplayPars(size=sizeTrack, color = "#00000080", type="line",lwd=1),trackOverlay=list(y.post.smooth))
	y.pre<-exprs(pSet[,t0])
	y.pre.smooth<-makeSmoothing(position(pSet), rowMeans(y.pre), dp = DisplayPars(color = "#80808090", type="line",lwd=4))
	y.pre.track<-makeGenericArray(as.matrix(y.pre), probeStart=position(pSet), dp = DisplayPars(size=sizeTrack, color = "#00000080", type="line",lwd=1),trackOverlay=list(y.pre.smooth))

	n.track<-4+(!minimalist)*4
	all.plot<-vector("list",n.track)
	all.plot[1]<-freq.track
	all.plot[2]<-anno.list[[1]]
	all.plot[3]<-y.track
	all.plot[4]<-anno.list[[1]]
	names(all.plot)[1:2]<-c("% Resp.","")
	names(all.plot)[3:4]<-c("Normalized\n Intensities\n log2(Post/Pre)","")
	
	if(!minimalist)
	{
	all.plot[5]<-y.post.track
	all.plot[6]<-anno.list[[1]]
	all.plot[7]<-y.pre.track
	all.plot[8]<-anno.list[[1]]
	names(all.plot)[5:6]<-c("Normalized\n log2(Post)","")
	names(all.plot)[7:8]<-c("Normalized\n log2(Pre)","")	
	}
	
	hotspots.highlight<-makeRectangleOverlay(start = start(hotspots), end = end(hotspots), region=c(1,length(anno.list)+n.track), dp = DisplayPars(color = "#FFFFBF",fill = "#FFFFBF", alpha = 0.2))

	pdPlot(c(anno.list, all.plot)
	, minBase=1, maxBase=857
	, labelRot=0
	, labelCex=.8
	, highlightRegions=hotspots.highlight)
}

.plot.clade<-function(pSet,freq,anno.list,minimalist=TRUE,colorTrack,colorAnno,sizeTrack,sizeAnno,hotspots)
{
	isLegend<-FALSE
	if(is.matrix(freq[[1]]))
	{
		isLegend<-TRUE
	}
	
	data.list<-lapply(1:7,function(i,x,y,dp,size,color){if(is.matrix(x[[i]])){color<-c("#00000090",paste(color[i],90,sep=""))}else{color<-color[i]};makeGenericArray(as.matrix(x[[i]]),probeStart=position(y[[i]]),dp=DisplayPars(size=size, color=color, type="line",lwd=4,cex=.2,axisCex=.7,isLegend=isLegend))},x=freq,y=pSet, size=sizeTrack, color = colorTrack,isLegend)

	all.plot<-vector("list",14)
	all.plot[seq(1,14,2)]<-data.list
	all.plot[seq(2,14,2)]<-anno.list[[1]]
	names(all.plot)<-rep("",14)
	names(all.plot)[seq(1,14,2)]<-c("M","A","B","C","D","CRF01","CRF02")

	pdPlot(c(anno.list
	   ,all.plot
	  )
	, minBase=1
	, maxBase=857
	, labelRot=0
	, labelCex=.6
	)
}