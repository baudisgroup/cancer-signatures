# Plot cytobands with high weightings from Autoencoder

rm(list=ls())
library(tidyverse)
library(karyoploteR)

breast_bands <- read_csv("/Users/bogao/DataFiles/new landscape/data/breast_bandsum.csv")


# the basic plot
pp <- getDefaultPlotParams(plot.type = 3)
title ='breast AE bands'

kp <- plotKaryotype(genome="hg38", chromosomes = unique(breast_bands$chro), plot.type=3, plot.params = pp,
                    labels.plotter = NULL, main=title)

kpDataBackground(kp)
kpDataBackground(kp, data.panel = 2)
kpAddChromosomeNames(kp, srt=45)
kpAxis(kp, r0=0.5, r1=1, tick.pos = c(0, 0.25, 0.5, 0.75, 1))
kpAxis(kp, r0=0.5, r1=0,  tick.pos = c(0.25, 0.5, 0.75, 1), labels = c(-0.25, -0.5, -0.75, -1))

kpAxis(kp, r0=0.5, r1=0, tick.pos = c(0, 0.25, 0.5, 0.75, 1), data.panel = 2)
kpAxis(kp, r0=0.5, r1=1,  tick.pos = c(0.25, 0.5, 0.75, 1), labels = c(-0.25, -0.5, -0.75, -1), data.panel = 2)


kp = kpBars(kp, chr=breast_bands$chro, x0=breast_bands$start, x1=breast_bands$end, y1=breast_bands$dup_scaled,
            r0=0.5, r1=1, col="#FFBD07", border=darker("#FFBD07"))

kp = kpBars(kp, chr=breast_bands$chro, x0=breast_bands$start, x1=breast_bands$end, y1=breast_bands$del_scaled,
            r0=0.5, r1=1, col="#00A6ED", border=darker("#00A6ED"), data.panel = 2)


# plot only top values
tops = c( top_n(breast_bands, 10, dup_scaled)$index, 
          top_n(breast_bands, 10, del_scaled)$index, 
          top_n(breast_bands, -10, dup_scaled)$index, 
          top_n(breast_bands, -10, del_scaled)$index)
breast_copy = filter(breast_bands, !index %in% tops)
breast_copy$dup_scaled = 0
breast_copy$del_scaled = 0
tops = filter(breast_bands, index %in% tops)
tops = bind_rows(breast_copy, tops)

kp <- plotKaryotype(genome="hg38", chromosomes = unique(breast_bands$chro), plot.type=3, plot.params = pp,
                    labels.plotter = NULL, main=title)

kpDataBackground(kp)
kpDataBackground(kp, data.panel = 2)
kpAddChromosomeNames(kp, srt=45)
kpAxis(kp, r0=0.5, r1=1, tick.pos = c(0, 0.25, 0.5, 0.75, 1))
kpAxis(kp, r0=0.5, r1=0,  tick.pos = c(0.25, 0.5, 0.75, 1), labels = c(-0.25, -0.5, -0.75, -1))

kpAxis(kp, r0=0.5, r1=0, tick.pos = c(0, 0.25, 0.5, 0.75, 1), data.panel = 2)
kpAxis(kp, r0=0.5, r1=1,  tick.pos = c(0.25, 0.5, 0.75, 1), labels = c(-0.25, -0.5, -0.75, -1), data.panel = 2)


kp = kpBars(kp, chr=tops$chro, x0=tops$start, x1=tops$end, y1=tops$dup_scaled,
            r0=0.5, r1=1, col="#3388FF", border="#3388FF")

kp = kpBars(kp, chr=tops$chro, x0=tops$start, x1=tops$end, y1=tops$del_scaled,
            r0=0.5, r1=1, col="#FE3231", border="#FE3231", data.panel = 2)  