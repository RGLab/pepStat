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

.plot.aggregate<-function(pSet,freq,anno.list,minimalist=TRUE,colorTrack,colorAnno,sizeTrack,sizeAnno,hotspots,ylim.freq=NULL)
{
    t0<-grepl("[Pp][Rr][Ee]",pData(pSet)$visit)
    t1<-grepl("[Pp][Oo][Ss][Tt]",pData(pSet)$visit)
	
	
	freq.track<-makeGenericArray(as.matrix(freq),probeStart=position(pSet),dp=DisplayPars(size=sizeTrack, color = paste(colorTrack[1:5],95,sep=""), type="line",lwd=4,cex=.2,axisCex=.7,legendCex=.7,isLegend=is.matrix(freq),ylim=ylim.freq))
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
	, labelCex=.6
	, highlightRegions=hotspots.highlight)
}

.plot.clade<-function(pSet,freq,anno.list,minimalist=TRUE,colorTrack,colorAnno,sizeTrack,sizeAnno,hotspots,ylim.freq=NULL)
{
	
	data.list<-lapply(1:7,function(i,x,y,dp,size,color,...){if(is.matrix(x[[i]])){color<-c("#00000090",paste(color[i],90,sep=""))}else{color<-color[i]};makeGenericArray(as.matrix(x[[i]]),probeStart=position(y[[i]]),dp=DisplayPars(legendCex=.7,size=size, color=color, type="line",lwd=4,cex=.2,axisCex=.7,isLegend=is.matrix(freq[[1]]),...))},x=freq,y=pSet, size=sizeTrack, color = colorTrack,ylim=ylim.freq)

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
	, labelCex=.7
	)
}

.plot.intensity<-function(y,y.smooth,position,anno.list,colorTrack,colorAnno,sizeTrack,sizeAnno,hotspots)
{

	if(is.list(y))
	{
		track.name<-paste(colnames(y[[1]]),c("(M)","(A)","(B)","(C)","(D)","(CRF01)","(CRF02)"))
		y<-lapply(y,function(x,y){colnames(x)<-y;x},"Normalized\n Intensities")
		y.smooth.track<-lapply(1:7,function(i,x,y,color){makeSmoothing(x[[i]],y[[i]], dp = DisplayPars(color = paste(color[i],"99",sep=""), type="line",lwd=6))},position,y.smooth,colorTrack)
		y.track<-lapply(1:7,function(i,x,y,z,...){makeGenericArray(as.matrix(x[[i]]), probeStart=y[[i]],trackOverlay=list(Smoothing=z[[i]]), ...)},y,position,y.smooth.track,dp = DisplayPars(size=sizeTrack, color=paste("#000000","75",sep=""), 
		type="dots",pointSize=.5,isLegend=T,legendCex=.6))
	
		n.track<-14
		all.plot<-vector("list",n.track)
		for(i in 1:7)
		{
		all.plot[i*2-1]<-y.track[[i]]
		all.plot[i*2]<-anno.list[[1]]
		names(all.plot)[(i*2-1):(i*2)]<-c(track.name[i],"")
		}
	}
	else
	{
	track.name<-colnames(y)
	colnames(y)<-"Normalized\n Intensities"
	y.smooth.track<-makeSmoothing(position, y.smooth, dp = DisplayPars(color = paste(colorTrack[2],"99",sep=""), type="line",lwd=6))
	y.track<-makeGenericArray(as.matrix(y), probeStart=position, dp = DisplayPars(size=sizeTrack, color=paste("#000000","75",sep=""), type="dots",pointSize=.5,isLegend=T,legendCex=.6),trackOverlay=list("Smoothing"=y.smooth.track))
	
	n.track<-2
	all.plot<-vector("list",n.track)
	all.plot[1]<-y.track
	all.plot[2]<-anno.list[[1]]
	names(all.plot)[1:2]<-c(track.name,"")
	}	
	hotspots.highlight<-NULL
	if(!is.null(hotspots))
	{
		hotspots.highlight<-makeRectangleOverlay(start = start(hotspots), end = end(hotspots), region=c(1,length(anno.list)+n.track), dp = DisplayPars(color = "#FFFFBF",fill = "#FFFFBF", alpha = 0.2))
	}

	pdPlot(c(anno.list,all.plot), minBase=1, maxBase=857, labelRot=0, labelCex=.7,highlightRegions=hotspots.highlight)
}

.plot.aggregate.freq<-function(pSet,freq,anno.list,colorTrack,colorAnno,sizeTrack,sizeAnno,ylim.freq=NULL)
{
    t0<-grepl("[Pp][Rr][Ee]",pData(pSet)$visit)
    t1<-grepl("[Pp][Oo][Ss][Tt]",pData(pSet)$visit)

	freq.track<-makeGenericArray(as.matrix(freq),probeStart=position(pSet),dp=DisplayPars(size=sizeTrack, color = paste(colorTrack[1:ncol(as.matrix(freq))],90,sep=""), type="line",lwd=6,cex=.2,axisCex=.7,legendCex=.7,isLegend=FALSE,ylim=ylim.freq))
	# freq.track<-makeGenericArray(as.matrix(freq),probeStart=position(pSet),dp=DisplayPars(size=sizeTrack, color = paste(colorTrack[1:2],95,sep=""), type="line",lwd=4,cex=.2,axisCex=.7,legendCex=.7,isLegend=FALSE,ylim=ylim.freq))

	n.track<-2
	all.plot<-vector("list",n.track)
	all.plot[1]<-freq.track
	names(all.plot)[1]<-"%\nresponders"
	all.plot[2]<-anno.list[[1]]
	names(all.plot)[2]<-""

	pdPlot(c(anno.list, all.plot)
	, minBase=1, maxBase=857
	, labelRot=0
	, labelCex=.9,
	highlightRegions=NULL)
}

.plot.clade.freq<-function(pSet,freq,anno.list,colorTrack,colorAnno,sizeTrack,sizeAnno,ylim.freq=NULL)
{
	
	data.list<-lapply(1:7,function(i,x,y,dp,size,color,...){color <- paste(color[1:ncol(as.matrix(x[[i]]))],90,sep="");makeGenericArray(as.matrix(x[[i]]),probeStart=position(y[[i]]),dp=DisplayPars(legendCex=.7,size=size, color=color, type="line",lwd=4,cex=.2,axisCex=.7,isLegend=F,...))},x=freq,y=pSet, size=sizeTrack, color = colorTrack,ylim=ylim.freq)

	all.plot<-vector("list",14)
	all.plot[seq(1,14,2)]<-data.list
	all.plot[seq(2,14,2)]<-anno.list[[1]]
	names(all.plot)<-rep("",14)
	names(all.plot)[seq(1,14,2)]<-c("M","A","B","C","D","CRF01","CRF02")

	pdPlot(c(anno.list,all.plot)
	, minBase=1
	, maxBase=857
	, labelRot=0
	, labelCex=.9
	)
}
