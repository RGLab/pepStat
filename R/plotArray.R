plotArray<-function(pSet,calls=NULL,array.ind=1,cName=NULL,paired=TRUE,pch=".",col=1,pch.calls="+",col.calls=2,show.inter=TRUE,show.union=FALSE,Anno=NULL,clade=NULL,...)
{

  if(!is.null(cName))
  {
    if(paired)
    {
      C<-exprs(pSet)[,grep(cName,sampleNames(pSet))]
    }
    else
    {
      
      C<-rowMeans(as.matrix(exprs(pSet)[,grep(cName,sampleNames(pSet))]))
    }
    lr<-exprs(pSet)[,-grep(cName,sampleNames(pSet))]-C
    sNames<-sampleNames(pSet)[-grep(cName,sampleNames(pSet))]
  }
  else
  {
    lr<-exprs(pSet)
    sNames<-sampleNames(pSet)
  }
  if(length(pSet@featurePosition)==0)
  {
    stop("Positions must be included in the pSet object!")
  }
  
  lr<-lr
  calls<-calls
  union<-apply(calls,1,any)
  inter<-apply(calls,1,all)
  
  if(paired & ncol(pData(pSet))>0)
  {
    trt<-pData(pSet)$trt[seq(1,ncol(pSet),2)]
  }
  else
  {
    trt<-pData(pSet)$trt
  }
  
  if(is.null(clade))
  {
    clade<-matrix(rep(TRUE,nrow(lr)),nrow(lr),1)
  }
  
  for(i in array.ind)
  {
    # ylim=c(min(lr[,i],-3),max(lr[,i]))
    ylim<-c(-3,9)
    main.title<-paste(sNames[i],", # of + calls:",sum(calls[,i]))
    if(ncol(pData(pSet))>0)
    {
      ## There is some meta data to be added
      main.title<-paste(main.title, "Trt:", trt[i])
    }
    
    plot(pSet@featurePosition[clade[,1]],lr[clade[,1],i],pch=pch,col=col[1],xlab="Position",ylab="Intensity",main=main.title,ylim=ylim,...)
    lines(smooth.spline(pSet@featurePosition[clade[,1]], lr[clade[,1],i]),lty=1,col=1, lwd=2)
    points(pSet@featurePosition[calls[,i] & clade[,1]],lr[,i][calls[,i] & clade[,1]], pch=pch.calls[1], col=col.calls[1])
    
    abline(h=0,col="grey",lwd=2)
    if(ncol(clade)>1)
    {
      for(j in 2:ncol(clade))
      {
        points(pSet@featurePosition[clade[,j]],lr[clade[,j],i],pch=pch,col=col[j],xlab="Position",ylab="Intensity",main=main.title,ylim=ylim,...)
        lines(smooth.spline(pSet@featurePosition[clade[,j]], lr[clade[,j],i]),lty=1,col=col[j],lwd=2)
        points(pSet@featurePosition[calls[,i] & clade[,j]],lr[,i][calls[,i] & clade[,j]], pch=pch.calls[j], col=col.calls[j])
      }
    }
    
    if(!is.null(Anno))
    {
      tmp<-apply(Anno,1,.AddAnno,min=0,max=ylim[2])
    }
    data(env)
    tmp<-apply(env[1:2,],1,.AddAnno2,min=ylim[1],max=ylim[1]+3)
    tmp<-apply(env[c(3:6,9),],1,.AddAnno3,min=ylim[1]+.1,max=ylim[1]+1.1)
    tmp<-apply(env[7:8,],1,.AddAnno3,min=ylim[1]+1.9,max=ylim[1]+2.9)
    tmp<-apply(env[10:11,],1,.AddAnno3,min=ylim[1]+.1,max=ylim[1]+1.1)
    
    if(show.union)
    {
      points(pSet@featurePosition[union],lr[,i][union],pch="v",col="blue",cex=1.2)
    }
    if(show.inter)
    {
      points(pSet@featurePosition[inter],lr[,i][inter],pch="^",col="blue",cex=1.2)
    }
    if(!all(clade))
    {
      legend(750,8,colnames(clade),lty=1,col=col, pch=pch.calls,lwd=2)
    }
  }
}

.AddAnno<-function(Anno,min,max)
{
  rect(Anno[2], min, Anno[3], max, density = NULL, angle = 45, col = NULL, border = NULL, lty = par("lty"), lwd = 1)
  # legend(as.double(Anno[2])-50,max,as.character(Anno[1]),bty="n",xjust=0,horiz=TRUE,text.width=1)
  text((as.double(Anno[2])+as.double(Anno[3]))/2, y = max-1, labels = as.character(Anno[1]), adj = NULL,pos = NULL, offset = 0.5, vfont = NULL,cex = .5, col = 4, font = NULL,srt=90)
  
}

.AddAnno2<-function(Anno,min,max)
{
  rect(Anno[2], min, Anno[3], max, density = 0, angle = -45, col = NULL, border = NULL, lty = par("lty"), lwd = 1)
  # legend(as.double(Anno[2])-50,max,as.character(Anno[1]),bty="n",xjust=0,horiz=TRUE,text.width=1)
  text(Anno[2], y = max-1, labels = as.character(Anno[1]), adj = NULL,pos = 4, offset = -.2, vfont = NULL,cex = .5, col = 4, font = NULL,srt=90)
  
}

.AddAnno3<-function(Anno,min,max)
{
  rect(Anno[2], min, Anno[3], max, density = 10, angle = 45, col = NULL, border = NULL, lty = par("lty"), lwd = 1)
  # legend(as.double(Anno[2])-50,max,as.character(Anno[1]),bty="n",xjust=0,horiz=TRUE,text.width=1)
  text(Anno[2], y = max-1, labels = as.character(Anno[1]), adj = NULL,pos = 4, offset = -0.2, vfont = NULL,cex = .5, col = 4, font = NULL,srt=90)
  
}


