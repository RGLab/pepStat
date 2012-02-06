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
    rownames(I)<-peptide(pSet)
	I	
}
