# Plot the bands signal after using AE and iNN weight for each disease or sub-group
library(tidyverse)
library(karyoploteR)
dd <- read_csv("/Users/bogao/DataFiles/new landscape/data/brain_sub.csv")

# the basic plot
pp <- getDefaultPlotParams(plot.type = 3)
pp$data1inmargin <- 1
pp$data2inmargin <- 1
pp$ideogramheight <- 5
pp$data1max=0.5
pp$data2max=0.5
pp$topmargin=50


title ='Brain 938* Glioma'

pdf('/Users/bogao/DataFiles/plots/newlandscape/brain938.pdf', width = 20, height=4)
kp <- plotKaryotype(genome="hg38", chromosomes = unique(dd$chro), plot.params = pp,  plot.type=3,
                    labels.plotter = NULL, main=title)

kpDataBackground(kp)
kpDataBackground(kp, data.panel = 2)
kpAddChromosomeNames(kp, srt=0, cex=0.6)

kp = kpBars(kp, chr=dd$chro, x0=dd$start, x1=dd$end, y1=dd$d938_dup,
            col="#3388FF", border=darker("#3388FF"))

kp = kpBars(kp, chr=dd$chro, x0=dd$start, x1=dd$end, y1=dd$d938_del,
            col="#FE3231", border=darker("#FE3231"),data.panel = 2)

dev.off()

plotBands = function(data, label, title, tofile=FALSE){
  # the basic plot
  pp <- getDefaultPlotParams(plot.type = 3)
  pp$data1inmargin <- 1
  pp$data2inmargin <- 1
  pp$ideogramheight <- 5
  pp$data1max=0.5
  pp$data2max=0.5
  pp$topmargin=50
  
  if (tofile){
  filepath = paste0('/Users/bogao/DataFiles/plots/newlandscape/brain', label, '.pdf')
  pdf(filepath, width = 20, height=4)
  }
  
  kp <- plotKaryotype(genome="hg38", chromosomes = unique(data$chro), plot.params = pp,  plot.type=3,
                      labels.plotter = NULL, main=title)
  
  kpDataBackground(kp)
  kpDataBackground(kp, data.panel = 2)
  kpAddChromosomeNames(kp, srt=0, cex=0.6)
  
  cname = paste0('d', label, '_dup')
  kp = kpBars(kp, chr=data$chro, x0=data$start, x1=data$end, y1=data[[cname]],
              col="#3388FF", border=darker("#3388FF"))
  cname = paste0('d', label, '_del')
  kp = kpBars(kp, chr=data$chro, x0=data$start, x1=data$end, y1=data[[cname]],
              col="#FE3231", border=darker("#FE3231"),data.panel = 2)
  
  if (tofile){
    dev.off()
  }
  
}

# plotBands(dd, '938', title='Brain 938* Glioma',tofile=TRUE)
# plotBands(dd, '944', title='Brain 944* Glioblastoma',tofile=TRUE)
# plotBands(dd, '940', title='Brain 940* Astroblastoma',tofile=TRUE)
# plotBands(dd, '945', title='Brain 945* Oligodendroglioma',tofile=TRUE)
# plotBands(dd, '939', title='Brain 939* Choroid plexus carcinoma',tofile=TRUE)
# plotBands(dd, '947', title='Brain 947* Medulloblastoma',tofile=TRUE)


