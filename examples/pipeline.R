## This example curated from the vignette -- please see vignette("pepStat")
## for more information
if (require("pepDat")) {

  ## Get example GPR files + associated mapping file
  dirToParse <- system.file("extdata/gpr_samples", package = "pepDat")
  mapFile <- system.file("extdata/mapping.csv", package = "pepDat")

  ## Make a peptide set
  pSet <- makePeptideSet(files = NULL, path = dirToParse,
                         mapping.file = mapFile, log=TRUE)

  ## Plot array images -- useful for quality control
  plotArrayImage(pSet, array.index = 1)
  plotArrayResiduals(pSet, array.index = 1, smooth = TRUE)

  ## Summarize peptides, using pep_hxb2 as the position database
  data(pep_hxb2)
  psSet <- summarizePeptides(pSet, summary = "mean", position = pep_hxb2)

  ## Normalize the peptide set
  pnSet <- normalizeArray(psSet)

  ## Smooth
  psmSet <- slidingMean(pnSet, width = 9)

  ## Make calls
  calls <- makeCalls(psmSet, freq = TRUE, group = "treatment",
                     cutoff = .1, method = "FDR", verbose = TRUE)

  ## Produce a summary of the results
  summary <- restab(psmSet, calls)

}
