makePeptideSetTest<-function(files=NULL, path=NULL, mapping.file=NULL, use.flags=FALSE,
                         rm.control.list=NULL,
                         norm.empty=FALSE, empty.control.list=c("empty","blank control"),
                         bgCorrect.method="normexp", log=TRUE, check.row.order=FALSE, verbose=FALSE)
{
  # There is some ambiguity with respect to what is Name and ID
  # ID -> peptide
  # Name -> annotation
  f <- function(x) as.numeric(x$Flags > -99)
  
  # before reading in files, check whether mapping.file is accessible,
  # to save user time in case they made a mistake.
  if (!is.null(mapping.file)){
    mapping.file<-.sanitize.mapping.file(mapping.file)
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
      error("No common features found across slides")
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
  if (log) {
    RG$E <- log2(RG$E)
    preproc(myDesc) <- c(preproc(myDesc),list(transformation="log", normalization="none"))
  } else {
    preproc(myDesc) <- c(preproc(myDesc),list(transformation="none", normalization="none"))
  }
  
  preproc(myDesc)$bgCorrect.method <- bgCorrect.method
  preproc(myDesc)$norm.empty <- norm.empty
  
  if (norm.empty) {
    #Check both name and ID for the control list
    index <- RG$genes$Name%in%empty.control.list | RG$genes$ID%in%empty.control.list
    if (verbose) {
      cat("** Background corrected using ", sum(index), " empty splots **\n")
    }
    mean.empty <- matrix(colMeans(as.matrix(RG$E[index,])), 
                         nrow=nrow(as.matrix(RG$E)), 
                         ncol=ncol(as.matrix(RG$E)), 
                         byrow=TRUE)
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
              featureRange = RangedData(IRanges(rep(0,nPep),rep(0,nPep)), featureID,
                                        peptide = featureSequence, block = layout$Block[ind.keep], row = layout$Row[ind.keep], column = layout$Column[ind.keep]), 
              exprs = as.matrix(RG$E-mean.empty)[ind.keep, ],
              experimentData = myDesc)
  
  # Make sure everything is stored as lower
  sampleNames(pSet)<-tolower(sampleNames(pSet))
  
  if(!is.null(mapping.file)){
    snamesIn <- sampleNames(pSet) %in% rownames(mapping.file)
    if(!all(snamesIn)){
      warning(paste(c("The following array samples were not found in mapping.file:",
                      sampleNames(pSet)[!snamesIn]), collapse = "\n"))
    }			
    pData(pSet) <- mapping.file[match(sampleNames(pSet), rownames(mapping.file)), ]
  }
  return(pSet)
}

.sanitize.mapping.file <- function(mapping.file)
{
  
  if (is.character(mapping.file)) {
    if (file.access(mapping.file, mode = 0) < 0)
      stop("mapping.file is not an accessible file path. Typo?")
    
    # check whether mapping file is a ".csv" file
    ext <- file_ext(mapping.file)
    if (ext != "csv")
      stop("Mapping file must be a .csv file")
    
    # ensure that mapping.file has a filename entry
    header <- scan(mapping.file, what = "character",
                           nlines = 1, sep = ",", quiet = TRUE)
    if( sum(c("filename","ptid","visit") %in% tolower(header)) != 3 )
      stop("mapping.file document header must include 3 mandatory columns: filename, ptid, visit")
    
    j <- match("filename", tolower(header))
    mapping.file <- read.csv(mapping.file, row.names=header[j])
  } else {
    if(!is.data.frame(mapping.file)){
      warning("Mapping.file object coerced to data frame")
      mapping.file <- as.data.frame(mapping.file)
    }
    header <- tolower(names(mapping.file))
    j <- match("filename", tolower(header))
    
    if (!is.na(j)) {
      row.names(mapping.file) <- mapping.file[,j]
      # drop the filename now that we are using it as row names
      mapping.file <- subset(mapping.file, select = -j)
    } else {
      message("filename entry not found in mapping.file, using rownames as file names")
    }
    if( sum(c("ptid", "visit") %in% tolower(colnames(mapping.file))) < 2 )
      stop("mapping.file object must include mandatory columns: ptid, visit")
  }
  
  rownames(mapping.file) <- tolower(rownames(mapping.file))
  colnames(mapping.file) <- tolower(colnames(mapping.file))
  
  mapping.file$ptid <- tolower(mapping.file$ptid)
  mapping.file$visit <- tolower(mapping.file$visit)		
  
  mapping.file
}