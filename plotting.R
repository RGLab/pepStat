library(pepStat)
library(Pviz)
load("big_restab.rda")

#opts
addNC <- TRUE
sel_clade <- LETTERS[1:3]

#Tracks
track_list <- list()
track_list <- c(track_list, ProteinAxisTrack(addNC = addNC))

if(length(sel_clade) > 0){
  for(clade in sel_clade){
    if(clade == sel_clade[length(sel_clade)]){
      legend <- TRUE
    } else{
      legend <- FALSE
    }
    track_list <- c(track_list, CladeTrack(rtl, clade, legend=legend, name = clade,
                                           type = "l"))
  }
}

#plot
plotTracks(track_list, from = 0, to = max(rtl$end))

