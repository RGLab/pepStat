readData<-function(files="",rm.control.list=c("empty","none","JPT-","Ig","Cy","landmark"))
{
  if(files=="")
  {
    files<-list.files(pattern=".xls")
  }
  else
  {
    files<-as.list(files)
  }
  message("Reading files: ", do.call(paste,as.list(files)))
  dataList<-lapply(files,function(x){data<-read.xls(x);data<-data[,which(!names(data)%in%c("Block","Column","Row"))];data})
  data<-do.call(cbind,dataList)
  data<-cbind(data$Name,data$Annotation,data[which(!names(data)%in%c("Name","Annotation"))])
  names(data)[1:2]<-c("Name","Annotation")
  
  #Need to check that the names and annotations are the same

  if(!is.null(rm.control.list))
  {
    ind.keep<-lapply(rm.control.list,function(x,Name){!grepl(x,Name)},as.character(data$Name))
    ind.keep<-do.call(cbind,ind.keep)
    ind.keep<-apply(ind.keep,1,all)
  }
  data<-data[ind.keep,]
  myDesc <- new("MIAME")
  preproc(myDesc)<-list(transformation="none", normalization="none")

  ### Setting the normalization parameters
  # Anno<-strsplit(as.character(data$Annotation),"_")
  # Anno<-do.call(rbind,Anno)
  # order<-order(as.integer(Anno[,1]),Anno[,2],Anno[,3])
  # newSet<-new('peptideSet',featureAnnotation=as.character(data$Annotation)[order], exprs=as.matrix(data[order,-c(1:2)]), featureSequence=as.character(data$Name)[order], experimentData=myDesc)
  newSet<-new('peptideSet',featureAnnotation=as.character(data$Annotation), exprs=as.matrix(data[,-c(1:2)]), featureSequence=as.character(data$Name), experimentData=myDesc)

}