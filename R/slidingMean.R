slidingMean <-function(peptideSet, width=5, verbose=FALSE, 
                       split.by.space=TRUE)
{
  if (class(peptideSet)!="peptideSet")
  {
    stop("peptideSet must be an object of class peptideSet")
  }
  if (preproc(peptideSet@experimentData)$transformation!="log" &
       preproc(peptideSet@experimentData)$transformation!="vsn") {
    stop("The probe measurements need to be log/vsn transformed!")
  }
  if (preproc(peptideSet@experimentData)$normalization=="none"){
    warning("You should probably normalize your data before using this function")
  }
  
  # Check whether user wishes to apply to multiple spaces 
  if (split.by.space & length(names(ranges(peptideSet))) > 1){
    oldrownames = rownames(ranges(peptideSet))
    s <- space(ranges(peptideSet))
    pSetList <- split(peptideSet, s)
    smoothedpSetList <- lapply(pSetList, function(set)
    {
      p <- position(set)
      o <- order(p)
      
      # reorder peptideSet by position order
      set <- set[o,]
      
      y <- exprs(set)
      p <- position(set)
      ny <- applySlidingMean(y, width, p)
      rownames(ny) <- rownames(y)[ro]
      exprs(set) <- ny
      rownames(exprs(set)) <- rownames(y)
      
      set
    })
    
    # collect the sorted RangedData into a single RangedData object,
    # put back into peptideSet
    rdata = do.call("rbind", lapply(smoothedpSetList, ranges))
    peptideSet@featureRange = rdata
    
    # collect the separate smoothed intensities into a single expression
    # matrix, gather rownames
    exprs(peptideSet) = do.call("rbind", lapply(smoothedpSetList, exprs))
    rownames(exprs(peptideSet)) = do.call("c", 
      lapply(smoothedpSetList, function(x) rownames(exprs(x))))
  }
  
  else {
    if (length(names(ranges(peptideSet))) > 1)
      warning("smoothing multiple spaces together in peptideSet object")
    # This could be made more efficient with multicore
    p <- position(peptideSet)
    o <- order(p)
    ro <- order(o)
    
    y <- exprs(peptideSet)[o,]
    p <- position(peptideSet)[o]
    ny <- applySlidingMean(y, width, p)
    exprs(peptideSet) <- ny[ro,]
  }
  
  if (verbose) {
    cat("** Finished processing ",nrow(peptideSet)," probes on ",ncol(peptideSet)," arrays **\n");
  }
  peptideSet
}


applySlidingMean <- function(y, width, position)
{
  yn <- sapply(position, function(p) {
    p.window <- abs(position - p) <= width/2
    colMeans(y[p.window,,drop = FALSE])
  })
  t(yn)
}



