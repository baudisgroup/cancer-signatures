rm(list=ls())
library(tidyverse)
library(karyoploteR)

plotAllGenes = function(title, rootpath='/Users/bogao/DataFiles/new landscape/files/', outpath=FALSE){
  # the basic plot
  pp <- getDefaultPlotParams(plot.type = 3)
  pp$data1inmargin <- 1
  pp$data2inmargin <- 1
  pp$ideogramheight <- 5
  pp$topmargin=50
  
  if (outpath != FALSE){
    pdf(outpath, width = 20, height=4)
  }
  
  amps = read_delim(paste0(rootpath, 'amp_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  dels = read_delim(paste0(rootpath, 'del_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  
  # amps = read_delim(paste0(rootpath, 'high_amp_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  # dels = read_delim(paste0(rootpath, 'high_del_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  
  # amp_census = read_delim(paste0(rootpath, 'amp_census.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  # del_census = read_delim(paste0(rootpath, 'del_census.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  
  kp <- plotKaryotype(genome="hg38", chromosomes = paste0('chr', seq(1:22)),  plot.params = pp,  plot.type=3,
                      labels.plotter = NULL, main=title)
  
  
  kpDataBackground(kp)
  kpDataBackground(kp, data.panel = 2)
  kpAddChromosomeNames(kp, srt=0, cex=0.6)
  
  ## all data use chro, individual disease uses chr
  ## all data use scaled, individual disease uses cnv_scaled
  
  # kp = kpBars(kp, chr=paste0('chr',amps$chro), x0=amps$start, x1=amps$end, y1=amps$scaled,
  #             col="#3388FF", border=darker("#3388FF"))
  kp = kpBars(kp, chr=paste0('chr',amps$chr), x0=amps$start, x1=amps$end, y1=amps$cnv_scaled,
              col="#3388FF", border=darker("#3388FF"))

  
  # kp = kpBars(kp, chr=paste0('chr',dels$chro), x0=dels$start, x1=dels$end, y1=-dels$scaled,
  #             col="#FE3231", border=darker("#FE3231"),data.panel = 2)
  kp = kpBars(kp, chr=paste0('chr',dels$chr), x0=dels$start, x1=dels$end, y1=-dels$cnv_scaled,
              col="#FE3231", border=darker("#FE3231"),data.panel = 2)


  
  if (outpath != FALSE){
    dev.off()
  }
  
}

# plotAllGenes('Overall Significant Alternations', outpath = '/Users/bogao/DataFiles/plots/newlandscape/analysis/all/all.pdf')
# plotAllGenes('Overall Significant Alternations', outpath = '/Users/bogao/DataFiles/plots/newlandscape/analysis/all_hi/all.pdf')
plotAllGenes('Overall Significant Alternations', 
             rootpath = '/Users/bogao/DataFiles/new landscape/files/new/',
             outpath = '/Users/bogao/DataFiles/plots/newlandscape/new/all.pdf')

##### plot all #####

outroot = '/Users/bogao/DataFiles/plots/newlandscape/analysis'
dirs = list.dirs('/Users/bogao/DataFiles/new landscape/files', recursive = FALSE)

for (subpath in dirs){
  organ = basename(subpath)
  diseases = list.dirs(subpath, recursive = FALSE)
  
  for (dis in diseases){
    disname = basename(dis)
    nameroot = paste(tolower(organ), disname, '', sep = '_')
    fpath = file.path(dis, nameroot)
    outpath = file.path(outroot, organ, disname)
    # print(file.path(outpath, 'all.pdf'))
    dir.create(outpath, showWarnings = FALSE,recursive = TRUE)
    plotAllGenes(title = paste(organ, disname, sep=': '),
                 rootpath = fpath,
                 outpath = file.path(outpath, 'all.pdf'))
  }
}

