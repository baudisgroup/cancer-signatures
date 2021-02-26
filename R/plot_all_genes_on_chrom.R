rm(list=ls())
library(tidyverse)
library(karyoploteR)

plotAllGenesOnChrom = function(title, chromosome, rootpath='/Users/bogao/DataFiles/new landscape/files/',  outpath=FALSE){
  # the basic plot
  pp <- getDefaultPlotParams(plot.type = 3)
  # pp$data1inmargin <- 1
  pp$data2inmargin <- 20
  pp$ideogramheight <- 40
  pp$topmargin=50
  
  if (outpath != FALSE){
    pdf(outpath, width = 20, height=4)
  }
  
  amps = read_delim(paste0(rootpath, 'amp_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  dels = read_delim(paste0(rootpath, 'del_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  amp_census = read_delim(paste0(rootpath, 'amp_census.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  del_census = read_delim(paste0(rootpath, 'del_census.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  
  kp <- plotKaryotype(genome="hg38", chromosomes = paste0('chr', chromosome),  plot.params = pp,  plot.type=3,
                      labels.plotter = NULL, main=title)
  
  kpAddCytobandLabels(kp, force.all = TRUE, srt=90,col="orange")
  kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=0.5,
                   minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")
  kpDataBackground(kp)
  kpDataBackground(kp, data.panel = 2)
  kpAddChromosomeNames(kp, srt=0, cex=0.8)
  
  ## all data use chro, individual disease uses chr
  ## all data use scaled, individual disease uses cnv_scaled
  
  # kp = kpBars(kp, chr=paste0('chr',amps$chro), x0=amps$start, x1=amps$end, y1=amps$scaled*0.8,
  #             col="#3388FF", border=darker("#3388FF"))
  # 
  # kp = kpPlotMarkers(kp, chr=paste0('chr',amp_census$chro), x=amp_census$start, labels=amp_census$symbol,
  #                   cex=0.5, y=0.85)
  # 
  # kp = kpBars(kp, chr=paste0('chr',dels$chro), x0=dels$start, x1=dels$end, y1=-dels$scaled*0.8,
  #             col="#FE3231", border=darker("#FE3231"),data.panel = 2)
  # 
  # kp = kpPlotMarkers(kp, chr=paste0('chr',del_census$chro), x=del_census$start, labels=del_census$symbol,data.panel = 2,
  #                    cex=0.5, y=0.85)
  
  kp = kpBars(kp, chr=paste0('chr',amps$chr), x0=amps$start, x1=amps$end, y1=amps$cnv_scaled*0.8,
              col="#3388FF", border=darker("#3388FF"))
  
  kp = kpPlotMarkers(kp, chr=paste0('chr',amp_census$chr), x=amp_census$start, labels=amp_census$symbol,
                     cex=0.5, y=0.85)
  
  kp = kpBars(kp, chr=paste0('chr',dels$chr), x0=dels$start, x1=dels$end, y1=-dels$cnv_scaled*0.8,
              col="#FE3231", border=darker("#FE3231"),data.panel = 2)
  
  kp = kpPlotMarkers(kp, chr=paste0('chr',del_census$chr), x=del_census$start, labels=del_census$symbol,data.panel = 2,
                     cex=0.5, y=0.85)
  
  
  if (outpath != FALSE){
    dev.off()
  }
  
}


# chroms = c(1,3,4,5,7,8,9,10,14,15,17,18,19,20,22)
chroms = c(4)
for (i in chroms){
  plotAllGenesOnChrom('Significant Alternations', i, outpath = paste0('/Users/bogao/DataFiles/plots/newlandscape/analysis/all/chr', i , '.pdf'))
}


# plotAllGenesOnChrom('Significant Alternations', 8)
# 
# plotAllGenesOnChrom('Significant Alternations', 9)
# plotAllGenesOnChrom('Significant Alternations', 10)
# plotAllGenesOnChrom('Significant Alternations', 18)

##### plot all #####

outroot = '/Users/bogao/DataFiles/plots/newlandscape/analysis'
dirs = list.dirs('/Users/bogao/DataFiles/new landscape/files', recursive = FALSE)
chroms = c(1:22)

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
    for (i in chroms){
      fname = paste0(i, '.pdf')
      x <- tryCatch(
        
        plotAllGenesOnChrom(title = paste(organ, disname, sep=': '), 
                            chromosome = i, 
                            rootpath = fpath,
                            outpath = file.path(outpath, fname))
        ,
        error = function(e){
          message(e)
          dev.off()
        }
      )
      
      

    }
  }
}

