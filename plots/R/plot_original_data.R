rm(list = ls())
library(tidyverse)
library(mongolite)
library(GenomicRanges)
library(karyoploteR)

###### original scripts #######

## Load data
m = mongo(collection = 'mecaned', db = 'Rebased')
data = m$find(
  query = '{}',
  fields = '{"_id":false, "normalized":false, "cytobands":false}')

meta = read_csv("/Users/bogao/DataFiles/new landscape/data/all_bands_meta.csv", col_types = "cccccdc")
d = filter(meta, site == 'Ovary' & morphology %in% c('80103', '84413', '84421')) %>% select('id')
ddd = left_join(d, data, by=c("id" = "sample_id"))

## all segments
dd = do.call(rbind, data$segments)

## filter normal data
thresh = 0.2
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

## plot
window = 1e5

pp <- getDefaultPlotParams(plot.type = 3)
pp$data1height <- 20
pp$data2height <- 20
pp$data1inmargin <- 1
pp$data2inmargin <- 1
pp$ideogramheight <- 2

kp <- plotKaryotype(genome="hg38", chromosomes = "autosomal", plot.params = pp,  plot.type=3, labels.plotter = NULL)

kpDataBackground(kp)
kpDataBackground(kp, data.panel = 2)
kpAddChromosomeNames(kp, srt=0, cex=0.6)

kp <- kpPlotDensity(kp, dups, col="#3388FF", border=darker("#3388FF"), window.size = window)
kp <- kpPlotDensity(kp, dels, data.panel = 2, col="#FE3231", border=darker("#FE3231"), window.size = window)

## plot one chromosome

kp <- plotKaryotype(genome="hg38", chromosomes = "chr1", plot.params = pp,  plot.type=3, labels.plotter = NULL)

kpDataBackground(kp)
kpDataBackground(kp, data.panel = 2)
kpAddChromosomeNames(kp, srt=0, cex=0.6)

kp <- kpPlotDensity(kp, dups, col="#3388FF", border=darker("#3388FF"), window.size = window)
kp <- kpPlotDensity(kp, dels, data.panel = 2, col="#FE3231", border=darker("#FE3231"), window.size = window)



################## Functions ################

## read all data from database
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

# plot the whole genome
plot_all = function(dups, dels, title='', window = 1e5, outpath=FALSE){

  pp <- getDefaultPlotParams(plot.type = 3)
  pp$data1height <- 20
  pp$data2height <- 20
  pp$data1inmargin <- 1
  pp$data2inmargin <- 1
  pp$ideogramheight <- 2
  
  # pp$data1inmargin <- 1
  # pp$data2inmargin <- 1
  # pp$ideogramheight <- 5
  # pp$topmargin=50
  
  if (outpath != FALSE){
    pdf(outpath, width = 20, height=4)
  }
  
  kp <- plotKaryotype(genome="hg38", chromosomes = "autosomal", plot.params = pp,  plot.type=3, labels.plotter = NULL, main = title)
  
  kpDataBackground(kp)
  kpDataBackground(kp, data.panel = 2)
  kpAddChromosomeNames(kp, srt=45, cex=1)
  
  kp <- kpPlotDensity(kp, dups, col="#3388FF", border=darker("#3388FF"), window.size = window)
  kp <- kpPlotDensity(kp, dels, data.panel = 2, col="#FE3231", border=darker("#FE3231"), window.size = window)
  
  if (outpath != FALSE){
    dev.off()
  }
}

# plot a chromosome
plot_chro = function(dups, dels, chr, title='', window = 1e5, outpath=FALSE){
  pp <- getDefaultPlotParams(plot.type = 3)
  # pp$data1height <- 20
  # pp$data2height <- 20
  pp$data2inmargin <- 20
  pp$ideogramheight <- 40
  pp$topmargin=50
  
  if (outpath != FALSE){
    pdf(outpath, width = 20, height=4)
  }
  
  kp <- plotKaryotype(genome="hg38", chromosomes = chr, plot.params = pp,  plot.type=3, labels.plotter = NULL, main = title)
  
  kpAddCytobandLabels(kp, force.all = TRUE, srt=90,col="orange")
  kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=0.5,
                   minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")
  kpDataBackground(kp)
  kpDataBackground(kp, data.panel = 2)
  kpAddChromosomeNames(kp, srt=0, cex=0.8)
  
  kp <- kpPlotDensity(kp, dups, col="#3388FF", border=darker("#3388FF"), window.size = window)
  kp <- kpPlotDensity(kp, dels, data.panel = 2, col="#FE3231", border=darker("#FE3231"), window.size = window)
  
  if (outpath != FALSE){
    dev.off()
  }
}

plot_disease = function(data, meta, site, morphology_list, title, outroot){
  dir.create(outroot, showWarnings = FALSE,recursive = TRUE)
  d = prep_data(data, meta, site, morphology_list)
  plot_all(d$dups, d$dels, title = title, outpath = paste0(outroot, 'all.pdf'))
  for (i in c(1:22)){
    plot_chro(d$dups, d$dels, chr=paste0('chr',i), title = title, outpath = paste0(outroot, i , '.pdf'))
  }
  
}

############## Test scripts ###########
plot_all(dups, dels)
plot_chro(dups, dels, 'chr2')

outroot = '/Users/bogao/DataFiles/plots/newlandscape/analysis/all/original/'
outroot = '/Users/bogao/DataFiles/plots/newlandscape/'
plot_all(dups, dels, title='All', outpath = paste0(outroot, 'all.pdf'))
plot_all(d$dups, d$dels, title='All', outpath = paste0(outroot, 'all.pdf'))
chroms = c(1,3,4,5,7,8,9,10,14,15,17,18,19,20,22)
for (i in chroms){
  plot_chro(dups, dels, chr = paste0('chr', i),  title='All', outpath = paste0(outroot, i , '.pdf'))
}

d = filter(data, topography == 'Ovary' & morphology %in% c('8010/3', '8441/3', '8442/1'))


t = prep_data(ddd)

plot_all(t$dups, t$dels)

data = load_data()
d = prep_data(data$data)

plot_all(d$dups, d$dels, title = "All", outpath = FALSE)
plot_chro(d$dups, d$dels, chr='chr8', title = "All")

d = prep_data(data$data, meta=data$meta, 'Cerebellum', c('80463'))
plot_all(d$dups, d$dels, title = "Non-small cell")
plot_chro(d$dups, d$dels, chr='chr5', title = "Non-small cell")


##### Run #####
diseases = list(
  list('Brain', 'Astrocytoma', c('94003', '94013')),
  list('Brain', 'Glioma', c('93803', '94403')),
  list('Brain', 'Mixed glioma', c('93823')),
  list('Brain', 'Oligodendroglioma', c('94503', '94513')),
  list('Brain', 'Primitive neuroectodermal tumor', c('94733')),
  
  list('Breast', 'Infiltrating duct carcinoma', c('85003')),
  list('Breast', 'Intraductal carcinoma', c('85002')),
  list('Breast', 'Lobular carcinoma', c('85203')),
  
  list('Colon', 'Adenocarcinoma', c('81403')),
  list('Colon', 'Adenocarcinoma intestinal type', c('81443')),
  list('Colon', 'Adenoma', c('81400')),
  list('Colon', 'Mucinous adenocarcinoma', c('84803')),
  
  list('Kidney', 'Clear cell adenocarcinoma', c('83103')),
  list('Kidney', 'Renal cell carcinoma', c('83123', '83173')),
  
  list('Liver', 'Hepatocellular carcinoma', c('81703')),
  
  list('Lung', 'Adenocarcinoma', c('81403', '82553')),
  list('Lung', 'Carcinoma', c('80103', '80123')),
  list('Lung', 'Non-small cell carcinoma', c('80463')),
  list('Lung', 'Small cell carcinoma', c('80413')),
  list('Lung', 'Squamous cell carcinoma', c('80703')),
  
  list('Ovary', 'Adenocarcinoma', c('81403', '83103', '83803')),
  list('Ovary', 'Carcinoma', c('80103', '84413', '84421')),
  list('Ovary', 'Mucinous cystadenoma', c('84700', '84800')),
  
  list('Prostate', 'Adenocarcinoma', c('81403')),
  
  list('Skin', 'Melanoma', c('87203', '87213', '87303')),
  
  list('Stomach', 'Adenocarcinoma', c('81403')),
  list('Stomach', 'Adenocarcinoma intestinal type', c('81443')),
  list('Stomach', 'Carcinoma diffuse type', c('81453')),
  list('Stomach', 'Gastrointestinal stromal sarcoma', c('89363')),
  list('Stomach', 'Tubular adenocarcinoma', c('82113'))
)


outroot = '/Users/bogao/DataFiles/plots/newlandscape/analysis'

for (dis in diseases){
  plot_disease(data$data, data$meta, site=dis[[1]], morphology_list = dis[[3]], title = paste(dis[[1]], dis[[2]], sep=': '),
             outroot = paste(outroot, dis[[1]], dis[[2]], '', sep='/'))
}



dis = list('Cerebellum', 'Medulloblastoma', c('94703', '94713', '94743'))
plot_disease(data$data, data$meta, site=dis[[1]], morphology_list = dis[[3]], title = paste(dis[[1]], dis[[2]], sep=': '),
             outroot = paste(outroot, dis[[1]], dis[[2]], '', sep='/'))