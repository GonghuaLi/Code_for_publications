options(warn = -1)
library(pheatmap)
library(rGPMM)
library(sybil)
library(ggsci)
library(ggplot2)
library(dplyr)
library(data.table)
library(EnhancedVolcano)

setwd('/home/ligh/projects/X/projects/Cell_aging/') #please use your own path

## construct pathway analysis
load('./data/cell_line_data/replication_senescence_data_v20231207.Rdata')
load('./data/cell_line_data/IRradiation_senescence_data_v20231207.Rdata')
load('./data/cell_line_data/ROS_senescence_data_v20231207.Rdata')
load('./data/cell_line_data/Oncogene_senescence_data_v20231207.Rdata')

SYBIL_SETTINGS("SOLVER", "cplexAPI")


sentypes = c('rep','ROS','IR','Onco')

prepross_GPMM <- function(flux3,clin){
  #flux3 = outlist$fluxes_by_Vmax_cor
  flux3.annote = Recon3.annote
  idxx = rowMeans(abs(flux3))> 1e-6 & rowSums(abs(flux3)>0) > 0.1*ncol(flux3)
  
  flux3 = flux3[idxx,]
  flux3.abs = fillgaps_rowMin(abs(flux3))
  
  idbb = CVs(flux3.abs) > 1e-3
  flux3.abs = flux3.abs[idbb,]  
  flux3 = flux3[idbb,]  
  #remove invalid mcmc model
  idyy = flux3.abs['biomass_reaction',] > 1e-9 
  flux3.abs = flux3.abs[,idyy]
  flux3 = flux3[,idyy]
  
  flux3.log2= as.matrix(log2(flux3.abs+ 1e-6))
  flux3.annote = flux3.annote[rownames(flux3),]
  clin.v = clin[colnames(flux3),]
  out = list()
  out[['flux3']] = flux3
  out[['flux3.log2']] = flux3.log2
  out[['flux3.annote']] = flux3.annote
  out[['clin.v']] = clin.v
  return(out)
}

for(i in 1:length(sentypes)){
  
  clin_name = paste0(sentypes[i],'Clin')
  expr_name = paste0(sentypes[i],'Expr')
  flux3vmax_name = paste0(sentypes[i],'Flux3vmax')
  xclin = get(clin_name)
  xExpr = get(expr_name)
  xflux = get(flux3vmax_name)
  xstudies = names(xclin)
  
  
  thenames = names(xflux)
  for(j in 1:length(thenames)){
    studyname = thenames[j]
    outdir1 = paste0('./results/',sentypes[i])
    if(!dir.exists(outdir1)){
      dir.create(outdir1)
    }
    rownames(xclin[[studyname]]) = colnames(xflux[[studyname]])
    gpmmresult = prepross_GPMM(xflux[[studyname]],xclin[[studyname]])
    outdir2 = paste0(outdir1,'/',studyname,'_vmax')
    print(outdir2)
    status = one_line_flux_analysis(gpmmresult$flux3,
                                    fluxCase = gpmmresult$flux3.log2[,gpmmresult$clin.v$ageType == 'sen'],
                                    fluxControl = gpmmresult$flux3.log2[,gpmmresult$clin.v$ageType == 'pro'],
                                    flux3.annote = gpmmresult$flux3.annote,
                                    outdir = outdir2)
  }
}