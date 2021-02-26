# Plot selected genes of selected diseases on selected chromosomes.

rm(list=ls())
library(tidyverse)
library(karyoploteR)
library(mongolite)


load_data = function(d='Rebased', c='mecaned'){
  m = mongo(collection = c, db = d)
  data = m$find(
    query = '{}',
    fields = '{"_id":false, "normalized":false, "cytobands":false}')
  
  meta = read_csv("/Users/bogao/DataFiles/new landscape/data/all_bands_meta.csv", col_types = "cccccdc")
  
  return(list("data"=data, "meta"=meta))
}

## filter data, seperate into dup/del and make genomeranges
prep_data = function(data, meta=FALSE, site=FALSE, morphology_list=FALSE, thresh=0.2){
  
  if (meta == FALSE){
    cohort_data = data
  } else {
    cohort = filter(meta, site == site & morphology %in% morphology_list) %>% select('id')
    cohort_data = left_join(cohort, data, by=c("id" = "sample_id"))
  }
  
  ## all segments
  dd = do.call(rbind, cohort_data$segments)
  
  ## filter normal data
  dd = dplyr::filter(dd, abs(value)>thresh & probes>3 & (end-start) >1)
  dd = mutate(dd, label = ifelse(value > 0, 'dup', 'del'))
  dd = dplyr::rename(dd, chromosome = chro)
  
  # reset levels
  dups = dplyr::filter(dd, label == 'dup') %>% makeGRangesFromDataFrame()
  newStyle = mapSeqlevels(seqlevels(dups), 'UCSC')
  dups = renameSeqlevels(dups, newStyle)
  
  dels = filter(dd, label == 'del') %>% makeGRangesFromDataFrame()
  newStyle = mapSeqlevels(seqlevels(dels), 'UCSC')
  dels = renameSeqlevels(dels, newStyle)
  
  return(list('dups'=dups, 'dels'=dels))
}

plotSelected = function(subpath, prefix, title, chros, amp_genes, del_genes, data, rootpath='/Users/bogao/DataFiles/new landscape/files/new/', 
                        window = 1e5, tofile=FALSE){
  # the basic plot
  pp <- getDefaultPlotParams(plot.type = 3)
  pp$data1inmargin <- 1
  pp$data2inmargin <- 1
  pp$ideogramheight <- 30
  pp$topmargin=50
  
  
  ## compact plot ##
  if (tofile){
    filepath = paste0(rootpath, subpath, 'gene_overlaps.pdf')
    pdf(filepath, width = 20, height=4)
  }
  

  
  
  
  # amps = read_delim(paste0(rootpath, subpath, prefix, 'amp_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  high_amps = read_delim(paste0(rootpath, subpath, prefix, 'high_amp_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  low_amps = read_delim(paste0(rootpath, subpath, prefix, 'low_amp_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)

  # dels = read_delim(paste0(rootpath, subpath, prefix, 'del_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  high_dels = read_delim(paste0(rootpath, subpath, prefix, 'high_del_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  low_dels = read_delim(paste0(rootpath, subpath, prefix, 'low_del_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)

  # kp <- plotKaryotype(genome="hg38", chromosomes = paste0('chr', seq(1:22)),  plot.params = pp,  plot.type=3,
  #                     labels.plotter = NULL, main=title)
  
  # chros = paste0('chr', c(4,5,6,7,8,9,10,11,12,14,16,17,18,19,20,22))
  kp <- plotKaryotype(genome="hg38", chromosomes = chros,  plot.params = pp,  plot.type=3,
                      labels.plotter = NULL, main=title)
  
  kpDataBackground(kp)
  kpDataBackground(kp, data.panel = 2)
  kpAddChromosomeNames(kp, srt=0, cex=1.6)
  
  
  kp <- kpPlotDensity(kp, data$dups, col="#EFF8FB", window.size = window)
  kp <- kpPlotDensity(kp, data$dels, data.panel = 2, col="#FBEFEF", window.size = window)
  
  
  # kp = kpBars(kp, chr=paste0('chr',amps$chr), x0=amps$start, x1=amps$end, y1=amps$cnv_scaled,
  #             col="#F1C40F", border=darker("#F1C40F"))
  kp = kpBars(kp, chr=paste0('chr',high_amps$chr), x0=high_amps$start, x1=high_amps$end, y1=high_amps$cnv_scaled,
              col="#3388FF", border=("#3388FF"))
  kp = kpBars(kp, chr=paste0('chr',low_amps$chr), x0=low_amps$start, x1=low_amps$end, y1=low_amps$cnv_scaled,
              col="#3388FF", border=darker("#3388FF"))
  kp = kpPlotMarkers(kp, chr=paste0('chr',amp_genes$chr), x=amp_genes$start, labels=amp_genes$symbol,
                     cex=1.3, data.panel = 2, y =0.4)
  
  # kp = kpBars(kp, chr=paste0('chr',dels$chr), x0=dels$start, x1=dels$end, y1=-dels$cnv_scaled,
  #             col="#F9EBEA", border=darker("#F9EBEA"),data.panel = 2)
  kp = kpBars(kp, chr=paste0('chr',high_dels$chr), x0=high_dels$start, x1=high_dels$end, y1=-high_dels$cnv_scaled,
              col="#FE3231", border=("#FE3231"),data.panel = 2)
  kp = kpBars(kp, chr=paste0('chr',low_dels$chr), x0=low_dels$start, x1=low_dels$end, y1=-low_dels$cnv_scaled,
              col="#FE3231", border=darker("#FE3231"),data.panel = 2)
  kp = kpPlotMarkers(kp, chr=paste0('chr',del_genes$chr), x=del_genes$start, labels=del_genes$symbol,
                     cex=1.3, y=0.4)

  
    
  if (tofile){
    dev.off()
  }
  
}



########################################
### Check the signature of subtypes ###
########################################
data = load_data()


chroms = paste0('chr', c(1,7,9,10,14,20))
genes = read_delim('/Users/bogao/DataFiles/Data/genome/protein_genes_biomart.tsv', "\t", escape_double = FALSE, trim_ws = TRUE,
                   col_names = c('id', 'description', 'chr', 'start', 'end', 'symbol'), skip=1)
## Glioma
# names = c('EGFR', 'MET', 'PDGFRA', 'MDM2', 'PIK3CA', 'CDK4', 'CDK6', 'CDKN2A', 'PTEN', 'RB1', 'KLF6', 'NF1', 'PPM1D', 'TP53', 'H3F3A',
#           'ST6GAL2',  'PIK3R1', 'CDKN2C', 'IDH1', 'ATRX')

amp_names = c('EGFR', 'MET', 'CDK6')
del_names = c('CDKN2A', 'PTEN', 'KLF6')
amp_genes = dplyr::filter(genes, symbol %in% amp_names)
del_genes = dplyr::filter(genes, symbol %in% del_names)
d = prep_data(data$data, data$meta, 'Brain', c('93803', '94403'))

plotSelected('Brain/Glioma/', 'brain_Glioma_', 'Brain Glioma: 9380/3, 9440/3',tofile = TRUE, 
                   chros = chroms, amp_genes=amp_genes, del_genes=del_genes, data=d)


## Melanoma
# names = c('NRAS', 'PTPRT','SALL4', 'CUX1', 'BRAF', 'SFRP4', 'RAC1', 'TRRAP', 'PPP6C', 'CDKN2A', 'GNAQ', 'XPA')
amp_names = c('NRAS', 'RAC1','TRRAP', 'PTPRT','SALL4', 'CUX1', 'BRAF', 'SFRP4')
del_names = c( 'PPP6C', 'CDKN2A', 'GNAQ', 'XPA')
amp_genes = dplyr::filter(genes, symbol %in% amp_names)
del_genes = dplyr::filter(genes, symbol %in% del_names)
d = prep_data(data$data, data$meta, 'Skin', c('87203', '87213', '87303'))

plotSelected('Skin/Melanoma/', 'skin_Melanoma_', 'Skin Melanoma: 8720/3',tofile = TRUE,
             chros = chroms, amp_genes=amp_genes, del_genes=del_genes, data=d)


## Medulloblastoma
# names = c('SUFU', 'ATM', 'KMT2C', 'PMS2', 'PTCH1','CTNNB1','DDX3X','SMARCA4','KMT2D','MYCN','BCOR','LDB1','GLI','OTX2','NOTCH','SNCAIP','SMO','PIK3CA')
amp_names = c('SMO',  'KMT2C', 'PMS2')
del_names = c('SUFU','PTCH1','LDB1','OTX2')
amp_genes = dplyr::filter(genes, symbol %in% amp_names)
del_genes = dplyr::filter(genes, symbol %in% del_names)
d = prep_data(data$data, data$meta, 'Cerebellum', c('94703', '94713', '94743'))

plotSelected('Cerebellum/Medulloblastoma/', 'cerebellum_Medulloblastoma_', 'Cerebellum Medulloblastoma: 9470/3, 9471/3, 9474/3',tofile = TRUE,
             chros = chroms, amp_genes=amp_genes, del_genes=del_genes, data=d)

