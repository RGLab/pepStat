#' Normalize tiling array data using sequence information
#'
#' This function is used to normalize the peptide microarray data using sequence
#' information.
#'
#' @param peptideSet A \code{peptideSet}. The expression data for the peptides as
#' well as annotations and ranges.
#' @param method A \code{character}. The normalization method to be used. Can be
#' "Zpep" or "ZpepQuad".
#' @param robust A \code{logical}. If TRUE, reweigthed least-squares estimates are
#'  computed.
#' @param centered A \code{logical}. If TRUE, recenter the data.
#'
#' @details
#' The available methods are "Zpep" and "ZpepQuad". These methods fit a linear model
#' using either linear or linear and quadratic terms (respectively), regressing
#' intensity on the peptides' five Z-scale scores. A peptide Z-scale score is
#' obtained by summing over the Z-scale values in Sandburg et al (1998) of the amino
#' acids the peptide comprises.
#'
#' Peptide Z-scale scores may be provided in the featureRange slot of peptideSet.
#' This slot is a \code{GRanges} object x, and the function will seek five columns
#' labelled z1 through z5 in values(x). If these are not found, the function attempts
#' to calculate Z-scales from sequence information found in peptide(peptideSet)
#'
#' If robust = TRUE the linear model is fit with t_4 distributed errors. The method
#' returns the residuals of each peptide intensity in the fitted linear model. If
#' centered = TRUE the fitted intercept term is added back to the residuals of the fit.
#'
#' @return A \code{peptideSet} object with updated normalized intensity values.
#'
#' @seealso \code{\link{summarizePeptides}}, \code{\link{makeCalls}}
#'
#' @author Raphael Gottardo, Gregory Imholte
#' @references Sandberg, M., Eriksson, L., Jonsson, J., Sjostrom, M., and Wold,
#' S. (1998). New chemical descriptors relevant for the design of biologically
#' active peptides. A multivariate characterization of 87 amino acids. Journal of
#'  Medicinal Chemistry 41, 2481-2491.
#'
#' @name normalizeArray
#' @rdname normalizeArray
#' @aliases NormalizeArray
#' @export
#' @example examples/pipeline.R
normalizeArray <- function(peptideSet, method = "ZpepQuad", robust = TRUE, centered = TRUE){
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
                               robust, centered)

  ### Set the normalized data as exprs
  exprs(peptideSet) <- ynew
  colnames(exprs(peptideSet)) <- sampleNames(peptideSet)

  ### Setting the normalization parameters
  if (robust){
    normalization <- paste(method, "robust", sep = " ")
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
  ind = grepl("z[1-9]", tolower(colnames(values(ranges(peptideSet)))))
  if(any(ind))
  {
    X <- sapply(which(ind), function(x) values(ranges(peptideSet))[[x]])
  }
  else
  {
    X <- .makeZpepMatrix(peptide(peptideSet))
  }

  if(method == "ZpepQuad")
    X <- cbind(X, X^2)

  as.matrix(X)
}
