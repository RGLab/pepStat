pat <- ProteinAxisTrack(addNC=TRUE)
ATprot <- ATrack(start=start(aprot), end=end(aprot), id=aprot$name, stacking="dense", fill=1:length(aprot), name="Proteins")
ATloop <- ATrack(start=start(aloop), end=end(aloop), id=aloop$name, stacking="dense", fill=length(aloop):1, name="Loops")
DT <- DTrack(start=start(ps), end=start(ps), data=V_calls[,2], type=c("l"), col=c("red", "yellow"))
PT <- ProbeTrack(sequence=peptide(ps), intensity=V_calls[,2], cex=0.2, probeStart=pepStat::position(ps), legend=TRUE, name="freq")
plotTracks(c(pat, ATprot, ATloop, PT), from=1, to=850, showFeatureId=TRUE)

apply(clade(ps), 2, function(x){ length(which(x == TRUE))})


load("tracks.rda")

#IN: ps, V_calls
clades <- colnames(clade(ps))
toPlot <- data.frame(position = upos)
for(curClade in clades){
  curPep <-  names(which(clade(ps)[, curClade] == TRUE))
  sps <- ps[curPep,]
  svc <- V_calls[curPep,]
  rownames(svc) <- rownames(sps) <- 1:nrow(sps)

  pos <- pepStat::position(sps)
  missingPos <- upos[!(upos %in% pos)]
  idcs <- sapply(missingPos, function(mp){
    idx <- max(which(pos - mp < 0))
  })
  missingvc <- ( svc[idcs,] + svc[idcs+1,] )/2

  data <- data.frame(rbind(svc, missingvc))
  colnames(data) <- paste(curClade, colnames(data), sep="_")
  data$position <- c(pos, missingPos)
  data <- data[ order(data$position), ]
  toPlot <- merge(toPlot, data, by="position")
}
head(toPlot)
dim(toPlot)

DT <- DTrack(start=toPlot$position, end=toPlot$position, data=toPlot[, 2:ncol(toPlot)], 
             groups=colnames(toPlot)[2:ncol(toPlot)], legend=TRUE)

DTP <- DTrack(start=toPlot$position, end=toPlot$position, data=toPlot[, grep("PLACEBO", colnames(toPlot), value=TRUE)],
              groups=grep("PLACEBO", colnames(toPlot), value=TRUE), legend=TRUE, name="Placebo group",
              cex=1.5, pch=1, type=c("p", "l"))
DTV <- DTrack(start=toPlot$position, end=toPlot$position, data=toPlot[, grep("VACCINE", colnames(toPlot), value=TRUE)],
              groups=grep("VACCINE", colnames(toPlot), value=TRUE), legend=TRUE, name="Vaccine group",
              cex=1.5, pch=1, type=c("p", "l"))

plotTracks(c(pat, ATprot, ATloop, DTP, DTV), from=1, to=850, showFeatureId=TRUE, background.title="darkgray")
plotTracks(c(pat, ATprot, ATloop, DTP, DTV), from=350, to=420, showFeatureId=TRUE)

dt <- data.table(V_calls)
dt[, position := pepStat::position(ps)]
dt <- dt[, lapply(.SD, mean), by="position"]
dt <- dt[order(position)]
dtl <- DTrack(start=dt$position, end=dt$position, data=dt[, colnames(dt)[2:ncol(dt)], with=FALSE],
              groups =  colnames(dt)[2:ncol(dt)], name="Frequence", legend=TRUE, type="l")

dtd <- DTrack()

loops <- getFeature(f, name=paste0("V", 1:5))
prots <- getFeature(f, category="protein")


dt <- DTrack(start=d2$position, end=d2$position, data=d2[, col])





data <- data.table(calls)
data <- data[, c("position", "clade") := list(position(peptideSet), ranges(peptideSet)$clade)]
#data <- data.table(melt(data, id.vars=c("clade", "position"), variable.name="group"))
ss <- strsplit(data$clade, ",")
nrep <- sapply(ss, length)
uclade <- unlist(ss)
data <- data[rep(1:nrow(data), nrep)]
data <- data[, clade := uclade]
data <- data[order(position)]

upos <- unique(data$position)
uclade <- unique(data$clade)
complete_data <- data.table(clade = rep(uclade, length(upos)), 
                            position = rep(upos, length(uclade)))




seq <- "asdsad-sadasaf--asdfasdf"
lenseq <- nchar(seq)
seq <- readLines("~/workspace/BioC/PEP.db/inst/extdata/alignments/Musclehxb2_7subtypes1L.fasta")[[2]]
f1 <- function(x){cumsum(sapply(unlist(strsplit(seq, "")), function(x){!grepl("-", x)}))}
f2 <- function(x){cumsum(sapply(unlist(strsplit(seq, "")), function(x){!grepl("-", x)}, USE.NAMES=FALSE ))}

f3 <- function(){
  pos <- c()
  cnt <- 0
  for(i in 1:nchar(seq)){
    if(substr(seq, i, i) == "-"){
    } else{
      cnt <- cnt+1
    }  
    pos <- c(pos, cnt)
  }
  return(pos)
}
