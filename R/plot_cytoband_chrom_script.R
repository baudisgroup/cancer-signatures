# Plot cytobands on a specific chromosome, which shows signal peaks in landscape plot.
# Testing script

rm(list=ls())
library(mongolite)
library(tidyverse)
library(karyoploteR)

m = mongo('bands_skin_normalized', 'Rebased')
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
title ='skin normalized'
# plot the chromsome bar and histograms
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes="chr11", plot.params = pp,labels.plotter = NULL, main=title)
kpDataBackground(kp)
kpDataBackground(kp, data.panel = 2)
kpAxis(kp)
kpAxis(kp,data.panel = 2)
kpAddChromosomeNames(kp)
kpAddCytobandLabels(kp, force.all=TRUE, srt=90,col="orange", cex=1)
kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_dup, col="#3388FF", border=darker("#3388FF"))
kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_del, col="#FE3231", border=darker("#FE3231"),data.panel = 2,r0=0, r1=-1)

# plot the original
m = mongo('bands_skin_original', 'Rebased')
bands = m$find()
bands = mutate(bands, ave_dup=total_dup/dup_length, ave_del=total_del/del_length)
gr_bands = makeGRangesFromDataFrame(dplyr::select(bands, chr,start,end,name,ave_dup,ave_del) %>%
                                      filter(!chr %in% c('chrX','chrY')), keep.extra.columns = TRUE)
title ='skin original'
kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_dup, col="#00E5EE", border=darker("#00E5EE"))
kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_del, col="#FFA54F", border=darker("#FFA54F"),data.panel = 2,r0=0, r1=-1)
