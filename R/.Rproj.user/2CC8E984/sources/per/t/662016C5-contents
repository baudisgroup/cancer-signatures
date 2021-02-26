# Plot selected genes of selected diseases on selected chromosomes.

rm(list=ls())
library(tidyverse)
library(karyoploteR)


names = c('EGFR', 'MET', 'PDGFRA', 'MDM2', 'PIK3CA', 'CDK4', 'CDK6', 'CDKN2A', 'PTEN', 'RB1', 'KLF6', 'NF1', 'PPM1D', 'TP53', 'H3F3A')
# getGenes = function(names){
genes = read_delim('/Users/bogao/DataFiles/Data/genome/protein_genes_biomart.tsv', "\t", escape_double = FALSE, trim_ws = TRUE )
s = dplyr::filter(genes, `HGNC symbol` %in% names)
names(s) = c('id', 'description', 'chr', 'start', 'end', 'symbol')
# }


plot_chro = function(subpath, prefix, title, chr, amp_genes, del_genes, rootpath='/Users/bogao/DataFiles/new landscape/files/', outpath=FALSE){
  pp <- getDefaultPlotParams(plot.type = 3)
  # pp$data1height <- 20
  # pp$data2height <- 20
  pp$data2inmargin <- 20
  pp$ideogramheight <- 40
  pp$topmargin=50
  
  if (outpath != FALSE){
    pdf(outpath, width = 20, height=4)
  }
  
  amps = read_delim(paste0(rootpath, subpath, prefix, 'amp_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  dels = read_delim(paste0(rootpath, subpath, prefix, 'del_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  
  kp <- plotKaryotype(genome="hg38", chromosomes = chr, plot.params = pp,  plot.type=3, labels.plotter = NULL, main = title)
  
  kpAddCytobandLabels(kp, force.all = TRUE, srt=90,col="orange")
  kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=0.5,
                   minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")
  kpDataBackground(kp)
  kpDataBackground(kp, data.panel = 2)
  kpAddChromosomeNames(kp, srt=0, cex=0.8)
  
  kp = kpBars(kp, chr=paste0('chr',amps$chr), x0=amps$start, x1=amps$end, y1=amps$cnv_scaled,
              col="#F1C40F", border=darker("#F1C40F"))
  kp = kpPlotMarkers(kp, chr=paste0('chr',amp_genes$chr), x=amp_genes$start, labels=amp_genes$symbol,
                     cex=1, data.panel = 2)
  
  
  kp = kpBars(kp, chr=paste0('chr',dels$chr), x0=dels$start, x1=dels$end, y1=-dels$cnv_scaled,
              col="#F9EBEA", border=darker("#F9EBEA"),data.panel = 2)
  kp = kpPlotMarkers(kp, chr=paste0('chr',del_genes$chr), x=del_genes$start, labels=del_genes$symbol,
                     cex=1)
  
  if (outpath != FALSE){
    dev.off()
  }
}




########################################
### Check the signature of subtypes ###
########################################




amp_names = c('EGFR', 'MET', 'PDGFRA', 'MDM2', 'PIK3CA', 'CDK4', 'CDK6')
del_names = c('CDKN2A', 'PTEN', 'RB1', 'KLF6', 'NF1', 'PPM1D', 'TP53', 'H3F3A')

# getGenes = function(names){
# }

genes = read_delim('/Users/bogao/DataFiles/Data/genome/protein_genes_biomart.tsv', "\t", escape_double = FALSE, trim_ws = TRUE,
                   col_names = c('id', 'description', 'chr', 'start', 'end', 'symbol'), skip=1)
amp_genes = dplyr::filter(genes, symbol %in% amp_names)
del_genes = dplyr::filter(genes, symbol %in% del_names)

plot_chro('Brain/Glioma/', 'brain_Glioma_', 'Brain Glioma: 9380/3, 9440/3', outpath = FALSE, 
           chr='chr10', amp_genes=amp_genes, del_genes=del_genes)



