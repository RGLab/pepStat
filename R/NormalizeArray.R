NormalizeArray = function(){
  .Deprecated("normalizeArray")
}

normalizeArray <- function(peptideSet, robust = TRUE, standard = FALSE,
                           method = "ZpepQuad", centered = TRUE)
{
  ### Sanity checks 
  if (class(peptideSet)!="peptideSet") {
    stop("peptideSet must be an object of class peptideSet (e.g. returned by makePeptideSet)!")
  }
  if (preproc(peptideSet@experimentData)$transformation != "log") {
    warning("The peptideSet measurements may not be log transformed!")
  }
  if (!(method %in% c("Zpep", "ZpepQuad")))
    stop("Invalid method argument")
  
  # Initialize Zpep
  Zpep <- getZpep(method, peptideSet)
  
  # Call fitting function
  ynew <- getWeightedEstimator(ymat = exprs(peptideSet), x = Zpep,
                               robust, centered, standard)
  
  ### Set the normalized data as exprs
  exprs(peptideSet) <- ynew
  colnames(exprs(peptideSet)) <- sampleNames(peptideSet)
  
  ### Setting the normalization parameters
  if (robust){
    normalization <- paste(method, "robust", sep = " ")
  }
  if (standard){
    normalization <- paste(method, "standardized", sep = " ")
  }
  preproc(peptideSet@experimentData)$normalization <- normalization
  peptideSet
}

getWeightedEstimator <- function(ymat, x, robust = TRUE, centered = TRUE,
                                 standard = FALSE) {
  x <- cbind(1, x)
  iter.max <- 10
  n <- nrow(x)
  out <- apply(ymat, 2, function(y) {
    # initial fit
    fit <- lm.fit(x = x, y = y)
    if (robust) {
      for(i in 1:iter.max) {
        RSS <- sum(fit$residuals^2)
        w <- (4 + 1)/(4 + (n/RSS)*fit$residuals^2)
        fit <- lm.wfit(x = x, y = y, w = w)
      }
    }
    if (standard) {
      o <- order(fit$fitted.values)
      ro <- order(o)
      
      n.bin <- 10
      n.probe <- length(y)
      n.ppb <- floor(n.probe/n.bin)
      rem <- n.probe - n.ppb*n.bin
      bin <- factor(c(rep(1:n.bin, each = n.ppb), rep(n.bin, rem)))
      
      split.resid <- split(fit$residuals[o], bin)
      scaled.resid <- unlist(lapply(split.resid, function(x) x/sd(x)))
      fit$residuals <- scaled.resid[ro]
    }
    if (!centered) {
      fit$residuals <- fit$residuals + fit$coefficients[1]
    }
    return(fit$residuals)
  })
  return(out)
}
