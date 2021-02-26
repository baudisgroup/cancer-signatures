# The function version to plot cytobands
# current version

rm(list=ls())
library(mongolite)
library(tidyverse)
library(karyoploteR)

# Plot the ideograph of data set and save a PDF to outpath
plotIdeos = function(db, collection, outpath){
  m = mongo(collection = collection, db = db)
  bands = m$find()
  
  bands = mutate(bands, ave_dup=total_dup/dup_length, ave_del=total_del/del_length)
  gr_bands = makeGRangesFromDataFrame(dplyr::select(bands, chr,start,end,name,ave_dup,ave_del)%>%filter(!chr %in% c('chrX','chrY')), keep.extra.columns = TRUE)
  
  pp <- getDefaultPlotParams(plot.type = 3)
  pp$data1height <- 20
  pp$data2height <- 20
  pp$data1inmargin <- 1
  pp$data2inmargin <- 1
  pp$ideogramheight <- 2
  
  # init pdf
  pdf(outpath, width = 20, height = 5)
  kp <- plotKaryotype(genome="hg38", plot.type=3, plot.params = pp,labels.plotter = NULL, main=collection)
  kpDataBackground(kp)
  kpDataBackground(kp, data.panel = 2)
  kpAxis(kp)
  kpAxis(kp,data.panel = 2)
  kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_dup, col="#3388FF", border=darker("#3388FF"))
  kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_del, col="#FE3231", border=darker("#FE3231"),data.panel = 2,r0=0, r1=-1)
  # save the file
  dev.off()
  
} 

# plot original and normalized for each dataset
projects = c('breast', 'skin', 'ovary')
dataform = c('original', 'normalized')
collection = apply(expand.grid('bands',projects, dataform), 1, paste, collapse="_")
for (c in collection){
  plotIdeos('Rebased', c, paste('/Users/bogao/DataFiles/plots/newlandscape/', c, '.pdf',sep ='' ))
}


# Combine the normalized and original data to one ideograph
plotCombinedIdeos = function(db, dataset, outpath){
  collection = paste(dataset,'normalized', sep='_')
  m = mongo(collection = collection, db = db)
  bands = m$find()
  
  bands = mutate(bands, ave_dup=total_dup/dup_length, ave_del=total_del/del_length)
  gr_bands = makeGRangesFromDataFrame(dplyr::select(bands, chr,start,end,name,ave_dup,ave_del)%>%filter(!chr %in% c('chrX','chrY')), keep.extra.columns = TRUE)
  
  pp <- getDefaultPlotParams(plot.type = 3)
  pp$data1height <- 20
  pp$data2height <- 20
  pp$data1inmargin <- 1
  pp$data2inmargin <- 1
  pp$ideogramheight <- 2
  
  pdf(outpath, width = 20, height = 5)
  
  # plot the normalized (bottom layer)
  kp <- plotKaryotype(genome="hg38", plot.type=3, plot.params = pp,labels.plotter = NULL, main=dataset)
  kpDataBackground(kp)
  kpDataBackground(kp, data.panel = 2)
  kpAxis(kp)
  kpAxis(kp,data.panel = 2)
  kpAddChromosomeNames(kp, srt=80)
  kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_dup, col="#3388FF", border=darker("#3388FF"))
  kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_del, col="#FE3231", border=darker("#FE3231"),data.panel = 2,r0=0, r1=-1)
  
  # prepare original data
  collection = paste(dataset,'original', sep='_')
  m = mongo(collection = collection, db = db)
  bands = m$find()
  
  bands = mutate(bands, ave_dup=total_dup/dup_length, ave_del=total_del/del_length)
  gr_bands = makeGRangesFromDataFrame(dplyr::select(bands, chr,start,end,name,ave_dup,ave_del)%>%filter(!chr %in% c('chrX','chrY')), keep.extra.columns = TRUE)
  
  # only need to plot data bars
  kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_dup, col="#00E5EE", border=darker("#00E5EE"))
  kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_del, col="#FFA54F", border=darker("#FFA54F"),data.panel = 2,r0=0, r1=-1)
  
  dev.off()
  
} 

# plot 3 projects.
collection = apply(expand.grid('bands',projects), 1, paste, collapse="_")
for (c in collection){
  plotCombinedIdeos('Rebased', c, paste('/Users/bogao/DataFiles/plots/newlandscape/', c, '.pdf',sep ='' ))
}
