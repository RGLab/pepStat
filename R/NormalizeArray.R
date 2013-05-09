NormalizeArray = function(peptideSet, robust = TRUE, standard = FALSE,
                          method = "ZpepQuad", centered = TRUE){
  .Deprecated("normalizeArray")
  normalizeArray(peptideSet, robust, standard, method, centered)
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


getZpep = function(method, peptideSet)
{
  ind = grepl("z[1-9]", tolower(colnames(ranges(peptideSet))))
  if(any(ind))
  {
    X <- sapply(which(ind), function(x) ranges(peptideSet)[[x]])
  }
  else
  {
    X <- makeZpepMatrix(peptide(peptideSet)) 
  }
  
  if(method == "ZpepQuad")
    X <- cbind(X, X^2)
  
  as.matrix(X)
}

computeZpep = function(AAstring, ztable)
{
  if(AAstring == c("empty"))
    return(rep(0, 5))
  t = unlist(strsplit(AAstring, split = ""))
  colSums(ztable[t,])
}

makeZpepMatrix = function(Sequence)
{
  Sequence = toupper(Sequence)
  let = unique(unlist(strsplit(Sequence, "")))
  AA = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
         "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  if(any(!(let %in% AA)))
    stop("Invalid sequences in Z-scale matrix construction")
  
  z = c(0.24, 0.84, 3.98, 3.11, -4.22,
        2.05, 2.47, -3.89, 2.29, -4.28,
        -2.85, 3.05, -1.66, 1.75, 3.52,
        2.39, 0.75, -2.59, -4.36, -2.54,
        -2.32, -1.67, 0.93, 0.26, 1.94,
        -4.06, 1.95, -1.73, 0.89, -1.3,
        -0.22, 1.62, 0.27, 0.5, 2.5,
        -1.07, -2.18, -2.64, 3.94, 2.44,
        0.6, 3.71, 1.93, -0.11, 1.06,
        0.36, 0.26, -1.71, -2.49, -1.49,
        0.47, 1.04, 1.84, -1.44, -3.5,
        1.15, -1.12, -1.54, 0.59, 0.43,
        -0.14, 0.18, -2.46, -3.04, 0.54,
        -0.82, 3.9, -0.84, 1.49, -0.72,
        1.94, -1.15, 0.7, -1.34, 1.99,
        -1.39, -1.46, -0.85, 3.44, 0.04,
        1.3, -2.65, 0.75, -0.25, -0.62,
        -0.38, 0.09, 0.26, 0.31, 0.84,
        -0.98, 1.61, 2, 0.66, -.17,
        0.67, -0.4, -0.02, -1.59, -1.47)
  dim(z) = c(20, 5)
  colnames(z) = paste0("z", 1:5)
  
  rownames(z) = AA
  
  Z = t(sapply(Sequence, computeZpep, ztable = z))
  colnames(Z) = paste0("z", 1:5)
  rownames(Z) = Sequence
  Z
}
