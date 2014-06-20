#' peptideSet constructor
#'
#' This function reads GenePix results (.gpr) files and creates a peptideSet object
#' gathering experiment information.
#'
#' @param files A \code{character} vector. If NULL, all files with a .gpr extension
#' in the specified path will be read.
#' @param path A \code{character} string. The directory where the .gpr files to
#' read are located.
#' @param mapping.file A \code{character} string or \code{data.frame}. A mapping file
#' that gives information for each sample. See details section below for a list of
#' required fields.
#' @param use.flags A \code{logical}. Should spots with flag value -99 or lower
#' be excluded?
#' @param rm.control.list A \code{character} vector. The name of the controls to
#' be excluded from the peptideSet.
#' @param empty.control.list A \code{character} vector. The name of the empty
#' controls. If non NULL, the intensity of these empty spots will be substracted
#' from remaining intensities.
#' @param bgCorrect.method A \code{character} string. The name of the method used
#' for background correction. This is passed to limma's backgroundCorrect method.
#' See details section below and ?backgroundCorrect for more information.
#' @param log A \code{logical}. If TRUE, intensities will be log2 transformed after
#' BG correction.
#' @param check.row.order A \code{logical}. Should slides be reduced to a common
#' set of peptides?
#' @param verbose A \code{logical}. Displays progress and additional information.
#'
#' @details
#' GenePix results files (.gpr) are read when found in either the files or path
#' arguments. By default, the foreground and background median intensities of the
#' "red" channels, "F635 Median" and "B635 Median", are read. The background
#' correction specified in bgCorrect.method is passed to the backgroundCorrect
#' method in the limma package.
#'
#' The mapping.file can be either a filename or a data.frame. In any case, it should
#' contain at least three columns labeled "filename", "ptid" and "visit". The
#' filenames given here should match those read from the path or files argument,
#' or be a subset of it. For each ptid (patient ID), the visit column should have at
#' least one "pre" and  one "post" sample. Any additional column will be kept into
#' the resulting \code{peptideSet} and can be used later on for groupping.
#'
#' If check.row.order = TRUE, the final set of probes is taken to be those with
#' IDs found in all arrays that were read.
#'
#' Row, Column and Block spatial array position for each probe are stored in the
#' \code{featureRanges} slot of the returned object.
#'
#'
#' @return A \code{peptideSet} object that contain the intensities, peptide
#' sequences and annotations available in the mapping file.
#'
#' @examples
#' # Read gpr files
#' library(pepDat)
#' mapFile <- system.file("extdata/mapping.csv", package = "pepDat")
#' dirToParse <- system.file("extdata/gpr_samples", package = "pepDat")
#' pSet <- makePeptideSet(files = NULL, path = dirToParse,
#'                        mapping.file = mapFile, log=TRUE)
#'
#' # Specify controls to be removed and empty controls
#' # to be used for background correction.
#' pSet <- makePeptideSet(files = NULL, path = dirToParse,
#'                        mapping.file = mapFile, log = TRUE,
#'                        rm.control.list = c("JPT-control", "Ig", "Cy3"),
#'                        empty.control.list= c("empty", "blank control"))
#'
#' @seealso \code{\link{peptideSet}}, \code{\link{read.maimages}},
#' \code{\link{backgroundCorrect}}
#'
#' @author Raphael Gottardo, Gregory Imholte
#'
#' @rdname makePeptideSet
#' @export
#' @importFrom tools file_path_sans_ext file_ext
#' @importFrom limma read.maimages backgroundCorrect
#' @importFrom Biobase preproc preproc<- assayData phenoData experimentData
#'   protocolData sampleNames sampleNames<- pData pData<- exprs exprs<- rowMedians
#' @importFrom IRanges IRanges space
#' @importFrom GenomicRanges GRanges
#'
makePeptideSet<-function(files=NULL, path=NULL, mapping.file=NULL, use.flags=FALSE,
                         rm.control.list=NULL, empty.control.list=NULL,
                         bgCorrect.method="normexp", log=TRUE, check.row.order=FALSE,
                         verbose=FALSE){
  # There is some ambiguity with respect to what is Name and ID
  # ID -> peptide
  # Name -> annotation
  f <- function(x) as.numeric(x$Flags > -99)

  # before reading in files, check whether mapping.file is accessible,
  # to save user time in case they made a mistake.
  if (!is.null(mapping.file)){
    mapping.file<-.sanitize_mapping_file2(mapping.file)
    files <- file_path_sans_ext(mapping.file$filename)
  }

  if (!check.row.order) { # Assume that the design is the same and don't check rows, order, etc.
    RG <- read.maimages(files=files,
                        source="genepix", path=path, ext="gpr",
                        columns=list(R="F635 Median",Rb="B635 Median"),
                        wt.fun=f,verbose=verbose)
  } else { # Used if the arrays don't exactly contain the same feature (e.g. the design has changed)
    files <- grep("gpr",list.files(path),value=TRUE)
    RG.list <- lapply(files, read.maimages, source="genepix",
                      path=path, columns=list(R="F635 Median",Rb="B635 Median"),
                      wt.fun=f, verbose=verbose)
    if(verbose){
      cat("Reordering all peptides\n")
    }

    RG <- RG.list[[1]]
    # Find the common target
    target.id <- Reduce(intersect,lapply(RG.list,function(x){x$genes$ID}))
    if(length(target.id)==0){
      stop("No common features found across slides")
    }

    # subset all
    RG.list <- lapply(RG.list,function(x, target){
      ind <- x$genes$ID%in%target;
      x$genes <- x$genes[ind,];
      x$Eb <- x$Eb[ind];
      x$E <- x$E[ind];
      x
    }, target.id)

    # order all
    RG.list <- lapply(RG.list, function(x, target){
      ind <- order(x$genes$ID);
      x$genes <- x$genes[ind,];
      x$Eb <- x$Eb[ind];
      x$E <- x$E[ind];
      x
    })
    RG$E <- do.call(cbind,lapply(RG.list,function(x){x$E}))
    RG$Eb <- do.call(cbind,lapply(RG.list,function(x){x$Eb}))
    RG$targets <- do.call(rbind,lapply(RG.list,function(x){x$targets}))
    # Make sure the sample names are consitent across objects
    colnames(RG$E) <- rownames(RG$targets)
    colnames(RG$Eb) <- rownames(RG$targets)
    RG$genes <- RG.list[[1]]$genes

  }

  offset <- 0.5
  if (bgCorrect.method=="half") offset <- .5 else offset <- 1
  RG <- try(backgroundCorrect(RG, method=bgCorrect.method, offset=offset, verbose=verbose))

  myDesc <- new("MIAME")

  ## Put NA instead of flags
  if (use.flags)
  {
    RG$E[RG$weights==0] <- NA
  }

  # Keep track of printer and source
  preproc(myDesc)$source<-RG$source
  preproc(myDesc)$printer<-RG$printer

  ## Fill in the details of the preprocessing

  ## Put NA instead of flags
  if (use.flags)
  {
    RG$E[RG$weights==0] <- NA
  }
  if (log) {
    RG$E <- log2(RG$E)
    preproc(myDesc) <- c(preproc(myDesc),list(transformation="log", normalization="none"))
  } else {
    preproc(myDesc) <- c(preproc(myDesc),list(transformation="none", normalization="none"))
  }

  if(!is.null(empty.control.list)){
    norm.empty <- TRUE
  } else{
    norm.empty <- FALSE
  }
  preproc(myDesc)$bgCorrect.method <- bgCorrect.method
  preproc(myDesc)$norm.empty <- norm.empty

  if (norm.empty) {
    #Check both name and ID for the control list
    index <- RG$genes$Name%in%empty.control.list | RG$genes$ID%in%empty.control.list
    if(sum(index)==0){
      warning("No empty controls matching the given list were found.")
      mean.empty <- rep(0, ncol(as.matrix(RG$E)))
    } else{
      if (verbose) {
        cat("** Background corrected using ", sum(index), " empty splots **\n")
      }
      mean.empty <- matrix(colMeans(as.matrix(RG$E[index,])),
                           nrow=nrow(as.matrix(RG$E)),
                           ncol=ncol(as.matrix(RG$E)),
                           byrow=TRUE)
    }
  } else {
    mean.empty <- rep(0, ncol(as.matrix(RG$E)))
  }

  # Keep the layout
  layout <- lapply(RG$genes[, c("Block","Row","Column")],as.factor)

  # Remove controls
  ind.keep <- rep(TRUE,nrow(RG$E))
  if (!is.null(rm.control.list)) {
    ind.keep <- lapply(rm.control.list,
                       function(x, Name, ID){
                         !grepl(x,Name) & !grepl(x,ID)},
                       as.character(RG$genes$Name),
                       as.character(RG$genes$ID))
    ind.keep <- do.call(cbind,ind.keep)
    ind.keep <- apply(ind.keep,1,all)
  }

  ## See note above about Name and ID
  # ID -> peptide
  # Name -> annotation
  featureSequence <- as.character(RG$genes$ID)[ind.keep]
  featureID <- as.character(RG$genes$Name)[ind.keep]
  nPep <- length(which(ind.keep))
  # Positions are set to zero before the information is provided in summarize pSet
  pSet <- new('peptideSet',
              featureRange = GRanges(seqnames = " ", strand = "*",
              ranges = IRanges(rep(0,nPep),rep(0,nPep)), featureID,
              peptide = featureSequence, block = layout$Block[ind.keep],
              row = layout$Row[ind.keep], column = layout$Column[ind.keep]),
              exprs = as.matrix(RG$E-mean.empty)[ind.keep, ],
              experimentData = myDesc)

  # Make sure everything is stored as lower
  sampleNames(pSet)<-tolower(sampleNames(pSet))

  if(!is.null(mapping.file)){
    pData(pSet) <- mapping.file[match(sampleNames(pSet), rownames(mapping.file)), ]
  }
  return(pSet)
}

.sanitize_mapping_file2 <- function(mapping.file){
  if(is.character(mapping.file)){
    map <- read.csv(mapping.file, header=TRUE)
  } else if(is.data.frame(mapping.file)){
    map <- mapping.file
  } else {
    stop("The mapping file should be a character vector or a data.frame")
  }
  colnames(map) <- tolower(colnames(map))
  if(!all(c("filename", "ptid", "visit") %in% colnames(map))){
    stop("The mapping file should include at least the 3 mandatory columns: 'filename', 'ptid' and 'visit'")
  }
  row.names(map) <- file_path_sans_ext(tolower(map$filename))
  return(map)
}
