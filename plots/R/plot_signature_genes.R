# Plot genes on genome for a disease data
rm(list=ls())
library(tidyverse)
library(karyoploteR)

plotSignatureGenes = function(subpath, prefix, title, rootpath='/Users/bogao/DataFiles/new landscape/files/new/', tofile=FALSE){
  # the basic plot
  pp <- getDefaultPlotParams(plot.type = 3)
  pp$data1inmargin <- 1
  pp$data2inmargin <- 1
  pp$ideogramheight <- 5
  pp$topmargin=50
  
  
  ## compact plot ##
  if (tofile){
    filepath = paste0(rootpath, subpath, 'sigGenes.pdf')
    pdf(filepath, width = 20, height=4)
  }
  
  high_amps = read_delim(paste0(rootpath, subpath, prefix, 'high_amp_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  low_amps = read_delim(paste0(rootpath, subpath, prefix, 'low_amp_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  amp_census = read_delim(paste0(rootpath, subpath, prefix, 'amp_census.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  
  high_dels = read_delim(paste0(rootpath, subpath, prefix, 'high_del_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  low_dels = read_delim(paste0(rootpath, subpath, prefix, 'low_del_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  del_census = read_delim(paste0(rootpath, subpath, prefix, 'del_census.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
    
  
  chros = paste0('chr', sort(unique(c(unique(high_amps$chr), unique(low_amps$chr),
                                      unique(high_dels$chr), unique(low_dels$chr)))))
  kp <- plotKaryotype(genome="hg38", chromosomes = chros,  plot.params = pp,  plot.type=3,
                      labels.plotter = NULL, main=title)
  
  kpDataBackground(kp)
  kpDataBackground(kp, data.panel = 2)
  kpAddChromosomeNames(kp, srt=0, cex=0.6)
  

  kp = kpBars(kp, chr=paste0('chr',high_amps$chr), x0=high_amps$start, x1=high_amps$end, y1=high_amps$cnv_scaled,
              col="#3388FF", border=("#3388FF"))
  kp = kpBars(kp, chr=paste0('chr',low_amps$chr), x0=low_amps$start, x1=low_amps$end, y1=low_amps$cnv_scaled,
              col="#3388FF", border=darker("#3388FF"))
  if (dim(amp_census)[1] >0){
    kp = kpPlotMarkers(kp, chr=paste0('chr',amp_census$chr), x=amp_census$start, labels=amp_census$symbol,
                    max.iter = 1000,cex=0.5, r0=0.25, y=0.25)
  }

  kp = kpBars(kp, chr=paste0('chr',high_dels$chr), x0=high_dels$start, x1=high_dels$end, y1=-high_dels$cnv_scaled,
              col="#FE3231", border=("#FE3231"),data.panel = 2)
  kp = kpBars(kp, chr=paste0('chr',low_dels$chr), x0=low_dels$start, x1=low_dels$end, y1=-low_dels$cnv_scaled,
              col="#FE3231", border=darker("#FE3231"),data.panel = 2)
  if (dim(del_census)[1] >0){
    kp = kpPlotMarkers(kp, chr=paste0('chr',del_census$chr), x=del_census$start, labels=del_census$symbol,data.panel = 2,
                     cex=0.5)
  }
  if (tofile){
    dev.off()
  }
  
  
  ## Full plot ##
  if (tofile){
    filepath = paste0(rootpath, subpath, 'sigGenes_full.pdf')
    pdf(filepath, width = 20, height=4)
  }
  
  high_amps = read_delim(paste0(rootpath, subpath, prefix, 'high_amp_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  low_amps = read_delim(paste0(rootpath, subpath, prefix, 'low_amp_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)

  high_dels = read_delim(paste0(rootpath, subpath, prefix, 'high_del_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)
  low_dels = read_delim(paste0(rootpath, subpath, prefix, 'low_del_genes.tsv'), "\t", escape_double = FALSE, trim_ws = TRUE)

  kp <- plotKaryotype(genome="hg38", chromosomes = "autosomal",  plot.params = pp,  plot.type=3,
                      labels.plotter = NULL, main=title)
  
  kpDataBackground(kp)
  kpDataBackground(kp, data.panel = 2)
  kpAddChromosomeNames(kp, srt=0, cex=0.6)
  
  kp = kpBars(kp, chr=paste0('chr',high_amps$chr), x0=high_amps$start, x1=high_amps$end, y1=high_amps$cnv_scaled,
              col="#3388FF", border=("#3388FF"))
  kp = kpBars(kp, chr=paste0('chr',low_amps$chr), x0=low_amps$start, x1=low_amps$end, y1=low_amps$cnv_scaled,
              col="#3388FF", border=darker("#3388FF"))

  kp = kpBars(kp, chr=paste0('chr',high_dels$chr), x0=high_dels$start, x1=high_dels$end, y1=-high_dels$cnv_scaled,
              col="#FE3231", border=("#FE3231"),data.panel = 2)
  kp = kpBars(kp, chr=paste0('chr',low_dels$chr), x0=low_dels$start, x1=low_dels$end, y1=-low_dels$cnv_scaled,
              col="#FE3231", border=darker("#FE3231"),data.panel = 2)
  
  
  if (tofile){
    dev.off()
  }
  
}

## Old code, merge code automaticaly
# plotSignatureGenes('Brain/', 'brain_', 'Brain All',tofile = TRUE)
# plotSignatureGenes('Brain/938/', 'brain_938_', 'Brain 938* Glioma',tofile = TRUE)
# plotSignatureGenes('Brain/939/', 'brain_939_', 'Brain 939* Choroid plexus carcinoma',tofile = TRUE)
# plotSignatureGenes('Brain/940/', 'brain_940_', 'Brain 940* Astroblastoma',tofile = TRUE)
# plotSignatureGenes('Brain/944/', 'brain_944_', 'Brain 944* Glioblastoma',tofile = TRUE)
# plotSignatureGenes('Brain/945/', 'brain_945_', 'Brain 945* Oligodendroglioma',tofile = TRUE)
# plotSignatureGenes('Brain/947/', 'brain_947_', 'Brain 947* Medulloblastoma',tofile = TRUE)
# 
# 
# plotSignatureGenes('Skin/', 'skin_', 'Skin All',tofile = TRUE)
# plotSignatureGenes('Skin/807/', 'skin_807_', 'Skin 807* Epidermoid/Squamous carcinoma',tofile = TRUE)
# plotSignatureGenes('Skin/872/', 'skin_872_', 'Skin 872* Melanoma',tofile = TRUE)
# plotSignatureGenes('Skin/873/', 'skin_873_', 'Skin 873* Amelanotic melanoma',tofile = TRUE)
# plotSignatureGenes('Skin/874/', 'skin_874_', 'Skin 874* Junctional Nevus',tofile = TRUE)
# plotSignatureGenes('Skin/877/', 'skin_877_', 'Skin 877* Mixed epithelioid and spindle cell melanoma',tofile = TRUE)
# plotSignatureGenes('Skin/883/', 'skin_883_', 'Skin 883* Undifferentiated Pleomorphic Sarcoma',tofile = TRUE)
# plotSignatureGenes('Skin/970/', 'skin_970_', 'Skin 970* Mycosis Fungoides',tofile = TRUE)
# plotSignatureGenes('Skin/C43/', 'skin_C43_', 'Skin C43* Melanoma',tofile = TRUE)


########################################
### Check the signature of subtypes ###
########################################

## First attempt, the final merged subgroup of OV
plotSignatureGenes('Ovary/', 'ovary_', 'Ovary All',tofile = TRUE)
plotSignatureGenes('Ovary/Carcinoma/', 'ovary_Carcinoma_', 'Ovary Carcinoma: 8010/3, 8441/3, 8442/1',tofile = TRUE)
plotSignatureGenes('Ovary/Adenocarcinoma/', 'ovary_Adenocarcinoma_', 'Ovary Adenocarcinoma: 8140/3, 8310/3, 8380/3',tofile = TRUE)
plotSignatureGenes('Ovary/Mucinous cystadenoma/', 'ovary_Mucinous cystadenoma_', 'Ovary Mucinous cystadenoma: 8470/0, 8480/0',tofile = TRUE)

## Brain: individual signatures
plotSignatureGenes('Brain/', 'brain_', 'Brain All',tofile = TRUE)
plotSignatureGenes('Brain/Primitive neuroectodermal tumor/', 'brain_Primitive neuroectodermal tumor_', 'Brain Primitive neuroectodermal tumor: 9473/3',tofile = TRUE)
plotSignatureGenes('Brain/Oligodendroglioma anaplastic/', 'brain_Oligodendroglioma anaplastic_', 'Brain Oligodendroglioma anaplastic: 9451/3',tofile = TRUE)
plotSignatureGenes('Brain/Oligodendroglioma/', 'brain_Oligodendroglioma_', 'Brain Oligodendroglioma: 9450/3',tofile = TRUE)
plotSignatureGenes('Brain/Glioblastoma/', 'brain_Glioblastoma_', 'Brain Glioblastoma: 9440/3',tofile = TRUE)
plotSignatureGenes('Brain/Astrocytoma anaplastic/', 'brain_Astrocytoma anaplastic_', 'Brain Astrocytoma anaplastic: 9401/3',tofile = TRUE)
plotSignatureGenes('Brain/Astrocytoma/', 'brain_Astrocytoma_', 'Brain Astrocytoma: 9400/3',tofile = TRUE)
plotSignatureGenes('Brain/Ependymoma/', 'brain_Ependymoma_', 'Brain Ependymoma: 9391/3',tofile = TRUE)
plotSignatureGenes('Brain/Mixed glioma/', 'brain_Mixed glioma_', 'Brain Mixed glioma: 9382/3',tofile = TRUE)
plotSignatureGenes('Brain/Glioma/', 'brain_Glioma_', 'Brain Glioma: 9380/3',tofile = TRUE)

## Brain: merged
plotSignatureGenes('Brain/', 'brain_', 'Brain All',tofile = TRUE)
plotSignatureGenes('Brain/Primitive neuroectodermal tumor/', 'brain_Primitive neuroectodermal tumor_', 'Brain Primitive neuroectodermal tumor: 9473/3',tofile = TRUE)
plotSignatureGenes('Brain/Oligodendroglioma/', 'brain_Oligodendroglioma_', 'Brain Oligodendroglioma: 9450/3, 9451/3',tofile = TRUE)
plotSignatureGenes('Brain/Astrocytoma/', 'brain_Astrocytoma_', 'Brain Astrocytoma: 9400/3, 9401/3',tofile = TRUE)
plotSignatureGenes('Brain/Ependymoma/', 'brain_Ependymoma_', 'Brain Ependymoma: 9391/3',tofile = TRUE)
plotSignatureGenes('Brain/Mixed glioma/', 'brain_Mixed glioma_', 'Brain Mixed glioma: 9382/3',tofile = TRUE)
plotSignatureGenes('Brain/Glioma/', 'brain_Glioma_', 'Brain Glioma: 9380/3, 9440/3',tofile = TRUE)

## Skin: individual signatures
plotSignatureGenes('Skin/', 'skin_', 'Skin All',tofile = TRUE)
plotSignatureGenes('Skin/Epidermoid carcinoma/', 'skin_Epidermoid carcinoma_', 'Skin Epidermoid carcinoma: 8070/3',tofile = TRUE)
plotSignatureGenes('Skin/Keratinizing/', 'skin_Keratinizing_', 'Skin Keratinizing: 8071/0',tofile = TRUE)
plotSignatureGenes('Skin/Melanoma/', 'skin_Melanoma_', 'Skin Melanoma: 8720/3',tofile = TRUE)
plotSignatureGenes('Skin/Nodular melanoma/', 'skin_Nodular melanoma_', 'Skin Nodular melanoma: 8721/3',tofile = TRUE)
plotSignatureGenes('Skin/Amelanotic melanoma/', 'skin_Amelanotic melanoma_', 'Skin Amelanotic melanoma: 8730/3',tofile = TRUE)
plotSignatureGenes('Skin/Bednar tumor/', 'skin_Bednar tumor_', 'Skin Bednar tumor: 8833/3',tofile = TRUE)
plotSignatureGenes('Skin/Pagetoid reticulosis/', 'skin_Pagetoid reticulosis_', 'Skin Pagetoid reticulosis: 9700/3',tofile = TRUE)

## Skin: merged
plotSignatureGenes('Skin/', 'skin_', 'Skin All',tofile = TRUE)
plotSignatureGenes('Skin/Epidermoid carcinoma/', 'skin_Epidermoid carcinoma_', 'Skin Epidermoid carcinoma: 8070/3',tofile = TRUE)
plotSignatureGenes('Skin/Keratinizing/', 'skin_Keratinizing_', 'Skin Keratinizing: 8071/0',tofile = TRUE)
plotSignatureGenes('Skin/Melanoma/', 'skin_Melanoma_', 'Skin Melanoma: 8720/3, 8721/3, 8730/3',tofile = TRUE)
plotSignatureGenes('Skin/Bednar tumor/', 'skin_Bednar tumor_', 'Skin Bednar tumor: 8833/3',tofile = TRUE)
plotSignatureGenes('Skin/Pagetoid reticulosis/', 'skin_Pagetoid reticulosis_', 'Skin Pagetoid reticulosis: 9700/3',tofile = TRUE)

## Breast: individual signatures
plotSignatureGenes('Breast/', 'breast_', 'Breast All',tofile = TRUE)
plotSignatureGenes('Breast/Pleomorphic carcinoma/', 'breast_Pleomorphic carcinoma_', 'Breast Pleomorphic carcinoma: 8022/3',tofile = TRUE)
plotSignatureGenes('Breast/Mucinous adenocarcinoma/', 'breast_Mucinous adenocarcinoma_', 'Breast Mucinous adenocarcinoma: 8480/3',tofile = TRUE)
plotSignatureGenes('Breast/Intraductal carcinoma/', 'breast_Intraductal carcinoma_', 'Breast Intraductal carcinoma: 8500/2',tofile = TRUE)
plotSignatureGenes('Breast/Infiltrating duct carcinoma/', 'breast_Infiltrating duct carcinoma_', 'Breast Infiltrating duct carcinoma: 8500/3',tofile = TRUE)
plotSignatureGenes('Breast/Intraductal micropapillary carcinoma/', 'breast_Intraductal micropapillary carcinoma_', 'Breast Intraductal micropapillary carcinoma: 8507/3',tofile = TRUE)
plotSignatureGenes('Breast/Lobular carcinoma/', 'breast_Lobular carcinoma_', 'Breast Lobular carcinoma: 8520/3',tofile = TRUE)
plotSignatureGenes('Breast/Infiltrating duct and lobular carcinoma/', 'breast_Infiltrating duct and lobular carcinoma_', 'Breast Infiltrating duct and lobular carcinoma: 8522/3',tofile = TRUE)
plotSignatureGenes('Breast/Infiltrating duct mixed/', 'breast_Infiltrating duct mixed_', 'Breast Infiltrating duct mixed: 8523/3',tofile = TRUE)
plotSignatureGenes('Breast/Inflammatory carcinoma/', 'breast_Inflammatory carcinoma_', 'Breast Inflammatory carcinoma: 8530/3',tofile = TRUE)
plotSignatureGenes('Breast/Metaplastic carcinoma/', 'breast_Metaplastic carcinoma_', 'Breast Metaplastic carcinoma: 8575/3',tofile = TRUE)

## Colon: individual signatures
plotSignatureGenes('Colon/', 'colon_', 'Colon All',tofile = TRUE)
plotSignatureGenes('Colon/Adenoma/', 'colon_Adenoma_', 'Colon Adenoma: 8140/0',tofile = TRUE)
plotSignatureGenes('Colon/Adenocarcinoma/', 'colon_Adenocarcinoma_', 'Colon Adenocarcinoma: 8140/3',tofile = TRUE)
plotSignatureGenes('Colon/Adenocarcinoma intestinal type/', 'colon_Adenocarcinoma intestinal type_', 'Colon Adenocarcinoma intestinal type: 8144/3',tofile = TRUE)
plotSignatureGenes('Colon/Mucinous adenocarcinoma/', 'colon_Mucinous adenocarcinoma_', 'Colon Mucinous adenocarcinoma: 8480/3',tofile = TRUE)

## Kidney: individual signatures
plotSignatureGenes('Kidney/', 'kidney_', 'Kidney All',tofile = TRUE)
plotSignatureGenes('Kidney/Oxyphilic adenoma/', 'kidney_Oxyphilic adenoma_', 'Kidney Oxyphilic adenoma 8290/0',tofile = TRUE)
plotSignatureGenes('Kidney/Clear cell adenocarcinoma/', 'kidney_Clear cell adenocarcinoma_', 'Kidney Clear cell adenocarcinoma 8310/3',tofile = TRUE)
plotSignatureGenes('Kidney/Renal cell carcinoma/', 'kidney_Renal cell carcinoma_', 'Kidney Renal cell carcinoma 8312/3',tofile = TRUE)
plotSignatureGenes('Kidney/Renal cell carcinoma chromophobe type/', 'kidney_Renal cell carcinoma chromophobe type_', 'Kidney Renal cell carcinoma chromophobe type 8317/3',tofile = TRUE)

## Liver: individual signatures
plotSignatureGenes('Liver/', 'liver_', 'Liver All',tofile = TRUE)
plotSignatureGenes('Liver/Hepatocellular carcinoma/', 'liver_Hepatocellular carcinoma_', 'Liver Hepatocellular carcinoma 8170/3',tofile = TRUE)

## Lung: individual signatures
plotSignatureGenes('Lung/', 'lung_', 'Lung All',tofile = TRUE)
plotSignatureGenes('Lung/Carcinoma/', 'lung_Carcinoma_', 'Lung Carcinoma 8010/3',tofile = TRUE)
plotSignatureGenes('Lung/Large cell carcinoma/', 'lung_Large cell carcinoma_', 'Lung Large cell carcinoma 8012/3',tofile = TRUE)
plotSignatureGenes('Lung/Small cell carcinoma/', 'lung_Small cell carcinoma_', 'Lung Small cell carcinoma 8041/3',tofile = TRUE)  
plotSignatureGenes('Lung/Non-small cell carcinoma/', 'lung_Non-small cell carcinoma_', 'Lung Non-small cell carcinoma 8046/3',tofile = TRUE)
plotSignatureGenes('Lung/Squamous cell carcinoma uncertain/', 'lung_Squamous cell carcinoma uncertain_', 'Lung Squamous cell carcinoma uncertain 8070/1',tofile = TRUE)
plotSignatureGenes('Lung/Squamous cell carcinoma/', 'lung_Squamous cell carcinoma_', 'Lung Squamous cell carcinoma 8070/3',tofile = TRUE)
plotSignatureGenes('Lung/Adenocarcinoma/', 'lung_Adenocarcinoma_', 'Lung Adenocarcinoma 8140/3',tofile = TRUE)
plotSignatureGenes('Lung/Bronchial adenoma carcinoid/', 'lung_Bronchial adenoma carcinoid_', 'Lung Bronchial adenoma carcinoid 8240/3',tofile = TRUE)
plotSignatureGenes('Lung/Bronchiolo-alveolar adenocarcinoma/', 'lung_Bronchiolo-alveolar adenocarcinoma_', 'Lung Bronchiolo-alveolar adenocarcinoma 8250/3',tofile = TRUE)
plotSignatureGenes('Lung/Bronchiolo-alveolar carcinoma non-mucinous/', 'lung_Bronchiolo-alveolar carcinoma non-mucinous_', 'Lung Bronchiolo-alveolar carcinoma non-mucinous 8252/3',tofile = TRUE)
plotSignatureGenes('Lung/Bronchiolo-alveolar carcinoma mixed mucinous and non-mucinous/', 'lung_Bronchiolo-alveolar carcinoma mixed mucinous and non-mucinous_', 'Lung Bronchiolo-alveolar carcinoma mixed mucinous and non-mucinous 8254/3',tofile = TRUE)
plotSignatureGenes('Lung/Adenocarcinoma with mixed subtypes/', 'lung_Adenocarcinoma with mixed subtypes_', 'Lung Adenocarcinoma with mixed subtypes 8255/3',tofile = TRUE)
### Papillary adenocarcinoma failed plotting
plotSignatureGenes('Lung/Papillary adenocarcinoma/', 'lung_Papillary adenocarcinoma_', 'Lung Papillary adenocarcinoma 8260/3',tofile = TRUE)
plotSignatureGenes('Lung/Mucinous adenocarcinoma/', 'lung_Mucinous adenocarcinoma_', 'Lung Mucinous adenocarcinoma 8480/3',tofile = TRUE)
plotSignatureGenes('Lung/Acinar cell carcinoma/', 'lung_Acinar cell carcinoma_', 'Lung Acinar cell carcinoma 8550/3',tofile = TRUE)
plotSignatureGenes('Lung/Adenosquamous carcinoma/', 'lung_Adenosquamous carcinoma_', 'Lung Adenosquamous carcinoma 8560/3',tofile = TRUE)

## Lung: merge
plotSignatureGenes('Lung/Adenocarcinoma/', 'lung_Adenocarcinoma_', 'Lung Adenocarcinoma 8140/3 8255/3',tofile = TRUE)
plotSignatureGenes('Lung/Carcinoma/', 'lung_Carcinoma_', 'Lung Carcinoma 8010/3, 8012/3',tofile = TRUE)

## Prostate: individual signatures
plotSignatureGenes('Prostate/', 'prostate_', 'Prostate All',tofile = TRUE)
plotSignatureGenes('Prostate/Carcinoma/', 'prostate_Carcinoma_', 'Prostate Carcinoma 8010/3',tofile = TRUE)
plotSignatureGenes('Prostate/Adenocarcinoma/', 'prostate_Adenocarcinoma_', 'Prostate Adenocarcinoma 8140/3',tofile = TRUE)

## Stomach: individual signatures
plotSignatureGenes('Stomach/', 'stomach_', 'Stomach All',tofile = TRUE)
plotSignatureGenes('Stomach/Carcinoma/', 'stomach_Carcinoma_', 'Stomach Carcinoma 8010/3',tofile = TRUE)
plotSignatureGenes('Stomach/Adenoma/', 'stomach_Adenoma_', 'Stomach Adenoma 8140/0',tofile = TRUE)
plotSignatureGenes('Stomach/Adenocarcinoma in situ/', 'stomach_Adenocarcinoma in situ_', 'Stomach Adenocarcinoma in situ 8140/2',tofile = TRUE)
plotSignatureGenes('Stomach/Adenocarcinoma/', 'stomach_Adenocarcinoma_', 'Stomach Adenocarcinoma 8140/3',tofile = TRUE)
plotSignatureGenes('Stomach/Adenocarcinoma intestinal type/', 'stomach_Adenocarcinoma intestinal type_', 'Stomach Adenocarcinoma intestinal type 8144/3',tofile = TRUE)
plotSignatureGenes('Stomach/Carcinoma diffuse type/', 'stomach_Carcinoma diffuse type_', 'Stomach Carcinoma diffuse type 8145/3',tofile = TRUE)
plotSignatureGenes('Stomach/Tubular adenocarcinoma/', 'stomach_Tubular adenocarcinoma_', 'Stomach Tubular adenocarcinoma 8211/3',tofile = TRUE)
plotSignatureGenes('Stomach/Mucinous adenocarcinoma/', 'stomach_Mucinous adenocarcinoma_', 'Stomach Mucinous adenocarcinoma 8480/3',tofile = TRUE)
plotSignatureGenes('Stomach/Signet ring cell carcinoma/', 'stomach_Signet ring cell carcinoma_', 'Stomach Signet ring cell carcinoma 8490/3',tofile = TRUE)
plotSignatureGenes('Stomach/Gastrointestinal stromal sarcoma/', 'stomach_Gastrointestinal stromal sarcoma_', 'Stomach Gastrointestinal stromal sarcoma 8936/3',tofile = TRUE)

### Cerebellum
plotSignatureGenes('Cerebellum/', 'cerebellum_', 'Cerebellum All',tofile = TRUE)
plotSignatureGenes('Cerebellum/Medulloblastoma/', 'cerebellum_Medulloblastoma_', 'Cerebellum Medulloblastoma 9470/3, 9471/3, 9474/3',tofile = TRUE)
# plotSignatureGenes('Cerebellum/Desmoplastic nodular medulloblastoma/', 'cerebellum_Desmoplastic nodular medulloblastoma_', 'Cerebellum Desmoplastic nodular medulloblastoma 9471/3',tofile = TRUE)
# plotSignatureGenes('Cerebellum/Large cell medulloblastoma/', 'cerebellum_Large cell medulloblastoma_', 'Cerebellum Large cell medulloblastoma 9474/3',tofile = TRUE)


