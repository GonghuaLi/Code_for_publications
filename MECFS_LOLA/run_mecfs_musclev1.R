setwd('/home/ligh/projects/X/projects/mecfs/')
library(rGPMM)
library(data.table)
library(GEOquery)

## muscle

prepross_GPMM <- function(outlist,clin){
  flux3 = outlist$fluxes_by_Vmax_cor
  flux3.annote = Recon3.annote
  idxx = rowMeans(abs(flux3))> 1e-6 & rowSums(abs(flux3)>0) > 0.1*ncol(flux3)
  
  flux3 = flux3[idxx,]
  flux3.abs = fillgaps_rowMin(abs(flux3))
  
  idbb = CVs(flux3.abs) > 1e-3
  flux3.abs = flux3.abs[idbb,]  
  flux3 = flux3[idbb,]  
  #remove invalid mcmc model
  idyy = flux3.abs['DM_atp_c_',] > 1e-9 
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


untar('./data/GSE245661/GSE245661_RAW.tar',exdir = './data/GSE245661/GSE245661')
files = list.files(path = './data/GSE245661/GSE245661',pattern = 'txt.gz',full.names = T)
tlist = list()
for(i in 1:length(files)){
  tmp = read.table(gzfile(files[i]))
  tmp  = tmp[-1,]
  tmp$V3 = as.numeric(tmp$V3)
  xx  = as.data.frame(t(tmp[,3]))
  colnames(xx) = tmp[,2]
  tlist[[i]] = xx
}
rawCounts = t(as.matrix(rbindlist(tlist,fill = T)))
ids = regmatches(files,regexpr('GSM\\d+',files))
colnames(rawCounts) = ids
fpkm = get_FPKM(rawCounts)
#rawCounts = cbind(data.frame(gene = rownames(rawCounts)),as.data.frame(rawCounts))
#

clin_org = getGEO(filename="./data/GSE245661/GSE245661_series_matrix.txt.gz")
clin = clin_org@phenoData@data
clin = clin[colnames(fpkm),]

changeSolver('cplex')

load('/home/ligh/projects/X/rGPMM/data/Recon3v2_GPMMmodel.Rdata')
load('/home/ligh/projects/X/rGPMM/data/reduceMod_3v2.Rdata')

idx = findRxnIDs(Recon3,'EX_nh4(e)')
Recon3$ub[idx] = 0
idx = findRxnIDs(Recon3,'EX_pyr(e)')
Recon3$ub[idx] = 0

GPMMdata$Recon3 = Recon3
GPMMdata$reduceMod = reduceMod


fluxout = rGPMM(expr = fpkm,expressionType = 'fpkm',standardization = F,tissueType = 'muscle')
gpmmresult = prepross_GPMM(fluxout,clin)

status = one_line_flux_analysis(gpmmresult$flux3,
                                fluxCase = gpmmresult$flux3.log2[,gpmmresult$clin.v$`subject status:ch1` != 'healthy volunteer (HV)'],
                                fluxControl = gpmmresult$flux3.log2[,gpmmresult$clin.v$`subject status:ch1` == 'healthy volunteer (HV)'],
                                flux3.annote = gpmmresult$flux3.annote,
                                outdir = './results/GSE245661v1/')


DEflux = file2frame('./results/GSE245661/DEflux.txt')
row.names(DEflux) = DEflux$ID

idx = DEflux$Pvalue < 0.05
DErxns = DEflux$ID[idx]
DEfc = DEflux$log2FC[idx]

# key metabolites
fpkm = as.data.frame(fpkm)
outdirkeymets = './results/GSE245661/keymetabolites.Rdata'
if(!file.exists(outdirkeymets)){
  keymetabolites = predict_key_metabolites(fpkm,DErxns,DEfc,numCores = 16)  #spend minutes
  save(list = c('keymetabolites'),file = outdirkeymets)
}
load('./results/GSE245661/keymetabolites.Rdata')
writetxt(x = keymetabolites$MetEffects,filename = './results/GSE245661/keymets.txt',row.names = T)

#key_genes
outdirkeygene = './results/GSE245661/keygenes.Rdata'
if(!file.exists(outdirkeygene)){
  keygenes = predict_key_genes(fpkm,DErxns,DEfc,numCores = 16)   #spend minutes
  save(list = c('keygenes'),file = outdirkeygene)
}

load('./results/GSE245661/keygenes.Rdata')
writetxt(x = keygenes$GeneEffects,filename = './results/GSE245661/keygenes.txt',row.names = T)



