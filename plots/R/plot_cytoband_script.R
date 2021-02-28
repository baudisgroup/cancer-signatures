# The general script to test and plot cytobands
# Can be used for test, not the final version.

rm(list=ls())
library(mongolite)
library(tidyverse)
library(karyoploteR)

# Easy reading data from mongodb to dataframe
m = mongo('bands_skin_original', 'Rebased')
bands = m$find()

# compute average values
bands = mutate(bands, ave_dup=total_dup/dup_length, ave_del=total_del/del_length)


##################################################################################
### Normalize all values, but it loses the meaning of values
### currently depreciated.

# vals = c(bands$ave_dup , -bands$ave_del)
# xmin = min(vals,na.rm = TRUE)
# xmax = max(vals,na.rm = TRUE)
# xdiff = xmax - xmin
# 
# bands = mutate(bands, ave_dup=(ave_dup-xmin)/xdiff, ave_del=(ave_del-xmin)/xdiff)
##################################################################################


# make a genomeRange 
gr_bands = makeGRangesFromDataFrame(dplyr::select(bands, chr,start,end,name,ave_dup,ave_del) %>%
                                      filter(!chr %in% c('chrX','chrY')), keep.extra.columns = TRUE)


# Plotting params, only needed once
pp <- getDefaultPlotParams(plot.type = 3)
pp$data1height <- 20
pp$data2height <- 20
pp$data1inmargin <- 1
pp$data2inmargin <- 1
pp$ideogramheight <- 2
title ='skin original'

# plot the chromsome bar and histograms
kp <- plotKaryotype(genome="hg38", plot.type=3, plot.params = pp,labels.plotter = NULL, main=title)
kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_dup, col="#00E5EE", border=darker("#00E5EE"))
kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_del, col="#FE3231", border=darker("#FE3231"),data.panel = 2,r0=0, r1=-1)


# Repeat the process for normalized data
m = mongo('bands_skin_normalized', 'Rebased')
bands = m$find()

bands = mutate(bands, ave_dup=total_dup/dup_length, ave_del=total_del/del_length)
gr_bands = makeGRangesFromDataFrame(dplyr::select(bands, chr,start,end,name,ave_dup,ave_del)%>%filter(!chr %in% c('chrX','chrY')), keep.extra.columns = TRUE)

title ='breast normalized'

# adding backgrounds and axis
kp <- plotKaryotype(genome="hg38", plot.type=3, plot.params = pp,labels.plotter = NULL, main=title)
kpDataBackground(kp)
kpDataBackground(kp, data.panel = 2)
kpAxis(kp)
kpAxis(kp,data.panel = 2)
kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_dup, col="#3388FF", border=darker("#3388FF"))
kp <- kpBars(kp, data=gr_bands, y0=0, y1=gr_bands$ave_del, col="#FE3231", border=darker("#FE3231"),data.panel = 2,r0=0, r1=-1)

