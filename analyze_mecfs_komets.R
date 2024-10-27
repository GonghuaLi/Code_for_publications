# plot ko metabolites

library(pheatmap)
library(rGPMM)
library(ggsci)
library(ggplot2)
library(dplyr)
library(data.table)
library(EnhancedVolcano)

get_status_from_ko <- function(KOeffects){
  ko = KOeffects
  tnames = names(ko)
  #del duplicated
  dupid = unique(tnames[duplicated(tnames)])
  if (length(dupid) > 0) {
    for(i in 1:length(dupid)){
      idx = which(tnames == dupid[i])
      tmpn = length(idx)
      if (any( ko[idx] > 0 ) & any(ko[idx] < 0)){
        ko[idx] = NA
      }else if(ko[idx[1]] > 0){
        ko[idx[-1]] = NA
      }else if(ko[idx[tmpn]] < 0){
        ko[idx[-tmpn]] = NA
      }else{
        ko[idx[-1]] = NA
      }
    }
    ko = ko[-which(is.na(ko))]
  }

  tsd = sd(ko)
  tpvalue = rep(1,length(ko))
  for(i in 1:length(tpvalue)){
    xx = ko[i]
    if(xx < 0){
      tpvalue[i] = pnorm(xx,0,tsd)
    }else{
      tpvalue[i] = 1-pnorm(xx,0,tsd)
    }
  }
  tfdr = p.adjust(tpvalue,method = 'fdr')
  out = data.frame(ID = names(ko),
                   ES = ko,
                   Pvalue = tpvalue,
                   FDR = tfdr)
  rownames(out) = out$ID
  return(out)
}

get_consist_mets <- function(ma){
  tmpout = ma
  tmpout$direct = sign(tmpout$ES)* (ma$FDR < 0.05 )
  tmpout$met_name = substr(tmpout$ID,1,nchar(tmpout$ID)-3)
  uname = unique(tmpout$met_name[tmpout$direct != 0])
  tmpout$direct.v = tmpout$direct
  for(i in 1:length(uname)){
    tmpids = which(tmpout$met_name == uname[i])
    tmpn = length(tmpids)
    if (tmpn < 2) {next;}
    tmps = tmpout$direct[tmpids]
    if(any(tmps > 0) & any(tmps < 0)){
      tmpout$direct.v[tmpids] = 0
    }else if(tmps[1] > 0){
      tmpout$direct.v[tmpids[-1]] = 0
    }else if(tmps[tmpn]< 0){
      tmpout$direct.v[tmpids[-tmpn]] = 0
    }else{
      tmpout$direct.v[tmpids[-1]] = 0
    }
  }
  return(tmpout)
}


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

load('./results/GSE245661/keymetabolites.Rdata')
tmp = keymetabolites$MetEffects
#KOmet[[i]] = get_status_from_ko(tmp)
tmp = keymetabolites$MetEffects
tmpkoxx = get_status_from_ko(tmp)
tmpkoxx = get_consist_mets(tmpkoxx)
#tmpkoxx = tmpkoxx[tmpkoxx$direct.v !=0,]
#rownames(tmpkoxx) = tmpkoxx$met_name
KOmet = tmpkoxx


tmpdata = data.frame(ID = rownames(KOmet),
                     log2FC = KOmet$ES,
                     Pvalue = as.vector(KOmet$FDR))
sum(tmpdata$Pvalue < 0.05 & tmpdata$log2FC > 0)
sum(tmpdata$Pvalue < 0.05 & tmpdata$log2FC < 0)

tmpdata$Pvalue[tmpdata$Pvalue < 1e-75] = 1e-75
tmpname = gsub('\\[\\w\\]','',tmpdata$ID)
cofactor = c('co2','o2','h2o','h','pi','nh4','nh3','no','ppi')
idx = !is.element(tmpname,cofactor)
tmpdata = tmpdata[idx,]

pdf('./Figures/KOmet_ES_new.pdf',width = 7,height = 7)
p = plot_DEflux(tmpdata,alpha = 0.3,num.showlab = 5,labSize = 7,fixpointsize = 2,
                xlab = bquote("Effective Score"),
                ylab = bquote(~-Log[10] ~ italic('FDR')),
                title = 'MECFS KOmets') + lghplot.addtheme(size = 24,sizex = 24,sizey = 24)
print(p)
dev.off()

colnames(tmpdata) = c('ID','ES score', 'FDR')
writetxt_forGPMM(tmpdata,filename = './results/GSE245661/keymetabolites_addfdr.txt')

koout = file2frame('./results/GSE245661/keymetabolites_addfdr.txt')
for(i in 1:nrow(koout)){
  idx = which(GPMMdata$Recon3$mets == koout$ID[i])
  koout$name[i] = GPMMdata$Recon3$metNames[idx]
}
writetxt_forGPMM(koout,'./results/GSE245661/keymetabolites_addfdr_A.txt')



