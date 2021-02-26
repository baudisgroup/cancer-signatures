# Test how to using karyoploteR and mongolite
# Test different plot types, see which suits my need.

rm(list=ls())
library(mongolite)
library(tidyverse)
library(karyoploteR)

# m = mongo('bands_breast_normalized', 'Rebased')
m = mongo('bands_skin_original', 'Rebased')
bands = m$find()

bands = mutate(bands, ave_dup=total_dup/dup_length, ave_del=total_del/del_length)
gr_bands = makeGRangesFromDataFrame(dplyr::select(bands, chr,start,end,name,ave_dup,ave_del)%>%filter(!chr %in% c('chrX','chrY')), keep.extra.columns = TRUE)



pp <- getDefaultPlotParams(plot.type = 3)
# pp$data1height <- 20
# pp$data2height <- 20
pp$data1inmargin <- 1
pp$data2inmargin <- 1
pp$ideogramheight <- 2
title ='breast normalized'

kp <- plotKaryotype(genome="hg38", plot.type=3, plot.params = pp, labels.plotter = NULL, main=title)
kpAddChromosomeNames(kp, srt=80)
kp <- kpLines(kp, data=gr_bands, y=gr_bands$ave_dup, col="#FFBD07", window.size = window)
kp <- kpLines(kp, data=gr_bands, y=gr_bands$ave_amp, data.panel = 2, col="#00A6ED", window.size = window)

kp <- kpPlotRibbon(kp, data=gr_bands, y0=0, y1=gr_bands$ave_dup, col="#FFBD07",r0=0, r1=-1)
kp <- kpPlotRibbon(kp, data=gr_bands, y0=0,y1=gr_bands$ave_amp, data.panel = 2, col="#00A6ED")

kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_dup, col="#00E5EE", border=darker("#00E5EE"),r0=0, r1=-1)

pp <- getDefaultPlotParams(plot.type = 3)
kp <- plotKaryotype(chromosomes="chr22", plot.type=2, plot.params = pp)
kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_dup, col="#00E5EE", border=darker("#00E5EE"),r0=0, r1=-1)

kp <- plotKaryotype(genome="hg38", plot.type=3, plot.params = pp,labels.plotter = NULL)
kpAxis(kp, tick.pos = c(0, 0.25, 0.5, 0.75, 1), labels=c(0,1,2,3,4))
kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_dup, col="#3388FF", border=darker("#3388FF"))
kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_del, col="#FE3231", border=darker("#FE3231"),data.panel = 2,r0=0, r1=-1)
