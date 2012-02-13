reportName<-tail(unlist(strsplit(getwd(),"/")),1)
path<-paste("~/Dropbox/Work/PeptideArray/Data","/",reportName,"/",sep="")
pathData<-paste(path,"gprFiles/",sep="")
width<-9
call.cutoff<-1.1
hotspot.cutoff<-20
imputeSaturated<-TRUE
minimalist<-FALSE