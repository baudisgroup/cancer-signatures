# Plot cytobands on a specific chromosome, which shows signal peaks in landscape plot.
# Functions version
# Current version

rm(list=ls())
library(mongolite)
library(tidyverse)
library(karyoploteR)

plotChrom = function(db, dataset, chromosome, outpath){
  
  collection = paste(dataset,'normalized', sep='_')
  m = mongo(collection, db)
  bands = m$find()
  
  # compute average values
  bands = mutate(bands, ave_dup=total_dup/dup_length, ave_del=total_del/del_length)
  
  gr_bands = makeGRangesFromDataFrame(dplyr::select(bands, chr,start,end,name,ave_dup,ave_del) %>%
                                        filter(!chr %in% c('chrX','chrY')), keep.extra.columns = TRUE)
  
  
  # Plotting params, only needed once
  pp <- getDefaultPlotParams(plot.type = 2)
  # pp$data1height <- 20
  # pp$data2height <- 20
  # pp$data1inmargin <- 1
  # pp$data2inmargin <- 1
  pp$ideogramheight <- 80
  
  pdf(outpath, width = 20, height = 5)
  
  # plot the chromsome bar and histograms
  kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes=chromosome, plot.params = pp,labels.plotter = NULL, main=collection)
  kpDataBackground(kp)
  kpDataBackground(kp, data.panel = 2)
  kpAxis(kp)
  kpAxis(kp,data.panel = 2)
  kpAddChromosomeNames(kp)
  kpAddCytobandLabels(kp, force.all=TRUE, srt=90,col="orange", cex=1)
  kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_dup, col="#3388FF", border=darker("#3388FF"))
  kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_del, col="#FE3231", border=darker("#FE3231"),data.panel = 2,r0=0, r1=-1)
  
  # plot the original
  collection = paste(dataset,'original', sep='_')
  m = mongo(collection, db)
  bands = m$find()
  bands = mutate(bands, ave_dup=total_dup/dup_length, ave_del=total_del/del_length)
  gr_bands = makeGRangesFromDataFrame(dplyr::select(bands, chr,start,end,name,ave_dup,ave_del) %>%
                                        filter(!chr %in% c('chrX','chrY')), keep.extra.columns = TRUE)
  kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_dup, col="#00E5EE", border=darker("#00E5EE"))
  kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_del, col="#FFA54F", border=darker("#FFA54F"),data.panel = 2,r0=0, r1=-1)
  dev.off()
}

# plots
plotChrom('Rebased','bands_breast', 'chr17', '/Users/bogao/DataFiles/plots/newlandscape/breast_17.pdf')
plotChrom('Rebased','bands_breast', 'chr8', '/Users/bogao/DataFiles/plots/newlandscape/breast_8.pdf')
plotChrom('Rebased','bands_breast', 'chr11', '/Users/bogao/DataFiles/plots/newlandscape/breast_11.pdf')
plotChrom('Rebased','bands_ovary', 'chr8', '/Users/bogao/DataFiles/plots/newlandscape/ovary_8.pdf')
plotChrom('Rebased','bands_ovary', 'chr11', '/Users/bogao/DataFiles/plots/newlandscape/ovary_11.pdf')
plotChrom('Rebased','bands_ovary', 'chr19', '/Users/bogao/DataFiles/plots/newlandscape/ovary_19.pdf')
plotChrom('Rebased','bands_skin', 'chr8', '/Users/bogao/DataFiles/plots/newlandscape/skin_8.pdf')
plotChrom('Rebased','bands_skin', 'chr11', '/Users/bogao/DataFiles/plots/newlandscape/skin_11.pdf')
plotChrom('Rebased','bands_skin', 'chr9', '/Users/bogao/DataFiles/plots/newlandscape/skin_9.pdf')
