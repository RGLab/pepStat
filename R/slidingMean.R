slidingMean <- function(peptideSet, width=9, verbose=FALSE){
  if (is(peptideSet,"peptideSet")){
    stop("peptideSet must be an object of class peptideSet")
  }
  if (preproc(peptideSet)$transformation!="log"){
    stop("The probe measurements need to be log/vsn transformed!")
  }
  if (preproc(peptideSet)$normalization=="none"){
    warning("You should probably normalize your data before using this function")
  }
  
  # Make a copy of the original ordering
  pep<-peptide(peptideSet)
  order<-order(position(peptideSet))
  
  # Reorder based on positions
  peptideSet<-peptideSet[order,]
  y <- exprs(peptideSet)
  
  if(all(position(peptideSet)==0)){
    stop("Your pSet object does not contain peptide position.\n Did you use `summarizePeptides' with a `position' argument?")
  }

  # This could be made more efficient with multicore
  ny<-apply(y, 2, slidingMeanVector, width, position(peptideSet))

  exprs(peptideSet)<-ny
  
  if(verbose)
  {
    cat("** Finished processing ",nrow(y)," peptides across ",ncol(y)," slides **\n");
  }
  # Reset the original ordering
  peptideSet[match(pep,peptide(peptideSet)),]
}

slidingMeanVector<-function(y, dMax, position)
{
  obj<-.C("MATScore",
  y,
  as.integer(length(y)),
  as.integer(1),
  as.integer(position),
  as.double(dMax),
  SScore = double(length(y)),
  as.integer(rep(1,length(y))),
  package = "pepStat")
  
  return(obj$SScore);
}
