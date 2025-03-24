setwd('/home/ligh/X/projects/MCMT/')
library(ggplot2)
library(reshape2)
library(Hmisc)
library(EnhancedVolcano)
library(pracma)
library(car)
library(ggrepel)
library(edgeR)
library(grid)
library(gridExtra)
library(Mfuzz)
library(M3C)
library(preprocessCore)
library(rlist)
library(ggsci)
library(scales)
library(data.table)
library(RColorBrewer)
library(readxl)
library(MetaDE)
#library(rGPMM)
library(RTN)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
source('./Subfunctions_Macaca_30tissue_aging.R')
set.seed(2025520)

library(msigdbr)
c2_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2")



updateRNA1 <- function(rnaData,headers,macacageneInfo){
  
  
  #filter with Rin > 6 and mapped counts > 10M
  idx = colSums(rnaData$rawCounts.whole) >= 10000000 & rnaData$rawCounts.whole.info$RIN >= 6
  rnaData$rawCounts.whole = rnaData$rawCounts.whole[,idx]
  rnaData$rawCounts.whole.info = rnaData$rawCounts.whole.info[idx,]
  
  outdata = rnaData
  
  tcounts = delete_dup_genes(rnaData$rawCounts.whole,headers)
  tcounts = tcounts[rownames(tcounts) != "",]
  
  #consider coding genes
  tmptype = macacageneInfo[rownames(tcounts),]$Type
  idx = tmptype == "protein_coding" | substr(tmptype,nchar(tmptype)-4, nchar(tmptype)) == "_gene"
  idx[is.na(idx)] = FALSE
  tcounts = tcounts[idx,]
  tcpm = get_CPM_counts(tcounts)
  tcpm.nonormal = get_CPM_nonormal(tcounts)
  
  batch = rnaData$rawCounts.whole.info$batch
  RINs = as.numeric(rnaData$rawCounts.whole.info$RIN)
  
  # 
  metadata <- data.frame(batch = batch, 
                         RIN = RINs,
                         log10counts = log10(colSums(rnaData$rawCounts.whole)),
                         age = rnaData$rawCounts.whole.info$age,
                         tissue = rnaData$rawCounts.whole.info$tissues)
  metadata$batch <- factor(metadata$batch)
  metadata$tissue <- factor(metadata$tissue)
  
  # Apply linear model to remove RINes
  #tcpm.batchad_pre = limma::removeBatchEffect(tcpm,covariates = metadata$RIN,
  #                                            design=model.matrix(~ age+tissue, data = metadata))
  #tcpm.batchad_pre = limma::removeBatchEffect(tcpm,covariates = as.matrix(metadata[,c("RIN","log10counts")]))
  tcpm.batchad_pre = limma::removeBatchEffect(tcpm,covariates = RINs)
  # Apply ComBat to remove batch effects
  tcpm.batchadj <- sva::ComBat(
    dat = tcpm.batchad_pre,
    mod = model.matrix(~ age+tissue, data = metadata),
    batch = batch,
    par.prior = TRUE,    # Use parametric adjustments
    prior.plots = FALSE  # Disable plotting of priors
  )
  
  
  #tissues
  tmp.mrna.tissues = list()
  tmp.mrna.tissues.info = list()
  tmp.mrna.tissues.org = list()
  tmp.mrna.tissues.org.info = list()
  tmp.tissue.names = unique(outdata$rawCounts.whole.info$tissues)
  for(i in 1:length(tmp.tissue.names)){
    thisname = tmp.tissue.names[i]
    vids = outdata$rawCounts.whole.info$tissues == thisname
    # get high confidance tissue data
    thiscouts = tcounts[,vids]
    tmpcpm.nonormal = tcpm.nonormal[,vids]
    xids = rowMeans(tmpcpm.nonormal) > 1 & rowSums(thiscouts < 5) < ncol(thiscouts)*0.2
    tmp.mrna.tissues[[thisname]] = tcpm.batchadj[xids,vids]
    tmp.mrna.tissues.info[[thisname]] = outdata$rawCounts.whole.info[vids,]
    tmp.mrna.tissues.org[[thisname]] = tcpm[xids,vids]
    tmp.mrna.tissues.org.info[[thisname]] = outdata$rawCounts.whole.info[vids,]
    
  }
  outdata$mrna.tissues = tmp.mrna.tissues
  outdata$mrna.tissues.org = tmp.mrna.tissues.org
  outdata$mrna.tissues.info = tmp.mrna.tissues.info
  return(outdata)
}


get_slid_mrnas_3_species_beta <- function(MetageneMouse,Metamrna,MetageneGETX){
  
  out = list()
  betacut.mouse = abs(sort(MetageneMouse$MetaBeta)[nrow(MetageneMouse)*0.1]) #0.1/18months
  betacut.macaca = abs(sort(Metamrna$MetaBeta)[nrow(Metamrna)*0.1]) #0.1/24years
  betacut.human = abs(sort(MetageneGETX$MetaBeta)[nrow(MetageneGETX)*0.1]) #0.1/60years
  
  #betacut.mouse = 0 #0.2/21months
  #betacut.macaca = 0 #0.2/27years
  #betacut.human = 0 #0.2/80years
  
  slids = (5:15)*100
  for(i in slids){
    tmpaa  = MetageneMouse[order(abs(MetageneMouse$MetaBeta),decreasing =T),]
    mouse_upmrna = tmpaa$ID[tmpaa$MetaBeta > betacut.mouse & tmpaa$manyNA == FALSE][1:i]
    mouse_downmrna = tmpaa$ID[tmpaa$MetaBeta < -betacut.mouse & tmpaa$manyNA == FALSE][1:i]
    
    
    tmpaa  = Metamrna[order(abs(Metamrna$MetaBeta),decreasing =T),]
    macaca_upmrna = tmpaa$ID[tmpaa$MetaBeta > betacut.macaca & tmpaa$manyNA == FALSE][1:i]
    macaca_downmrna = tmpaa$ID[tmpaa$MetaBeta < -betacut.macaca & tmpaa$manyNA == FALSE][1:i]
    
    
    tmpaa  = MetageneGETX[order(abs(MetageneGETX$MetaBeta),decreasing =T),]
    human_upmrna = tmpaa$ID[tmpaa$MetaBeta > betacut.human & tmpaa$manyNA == FALSE][1:i]
    human_downmrna = tmpaa$ID[tmpaa$MetaBeta < -betacut.human & tmpaa$manyNA == FALSE][1:i]
    out[[paste0('TOP',i)]] = list(mouse_upmrna = mouse_upmrna,
                                  mouse_downmrna = mouse_downmrna,
                                  macaca_upmrna = macaca_upmrna,
                                  macaca_downmrna = macaca_downmrna,
                                  human_upmrna = human_upmrna,
                                  human_downmrna = human_downmrna)
  }
  return(out)
  
}



get_GO_slid_mrnas_1 <- function(slid_mrnas,TERM2GENE){
  slids = (5:15)*100
  
  out = list()
  
  for(i in slids){
    xid = paste0('TOP',i)
    print(xid)
    tmpmrna = slid_mrnas[[xid]]
    mouse_upmrna  = GOenrichment.C2(tmpmrna$mouse_upmrna,TERM2GENE = TERM2GENE)
    mouse_downmrna  = GOenrichment.C2(tmpmrna$mouse_downmrna,TERM2GENE = TERM2GENE)
    macaca_upmrna  = GOenrichment.C2(tmpmrna$macaca_upmrna,TERM2GENE = TERM2GENE)
    macaca_downmrna  = GOenrichment.C2(tmpmrna$macaca_downmrna,TERM2GENE = TERM2GENE)
    human_upmrna  = GOenrichment.C2(tmpmrna$human_upmrna,TERM2GENE = TERM2GENE)
    human_downmrna  = GOenrichment.C2(tmpmrna$human_downmrna,TERM2GENE = TERM2GENE)
    out[[xid]] = list(mouse_upmrna = mouse_upmrna,
                                  mouse_downmrna = mouse_downmrna,
                                  macaca_upmrna = macaca_upmrna,
                                  macaca_downmrna = macaca_downmrna,
                                  human_upmrna = human_upmrna,
                                  human_downmrna = human_downmrna)
  }
  
  return(out)
}


plot_GOenrich_ratio_1 <- function(slid_mrnas,slid_GO,outfile){
  
  #### overlap genes
  slids = (5:15)*100
  k = 0;
  out.rho.up = matrix(0,3,length(slids))
  colnames(out.rho.up)  =  paste0('TOP',slids)
  rownames(out.rho.up)  = c('C2_macaca_vs_mouse','C2_macaca_vs_human','C2_human_vs_mouse')
  
  out.rho.down = out.rho.up
  pdf(outfile,width = 8)
  for(i in slids){
    xid = paste0('TOP',i)
    tmpmrna= slid_mrnas[[xid]]
    tmpGO = slid_GO[[xid]]
    
    p = ggvenn::ggvenn(data = list(mrna.up.human = tmpmrna$human_upmrna,
                                   mrna.up.macaca = tmpmrna$macaca_upmrna,
                                   mrna.up.mouse = tmpmrna$mouse_upmrna),
                       fill_color = c("#E64B35FF", "#3C5488FF", "green"),
                       show_percentage = F,stroke_size = 0.5,
                       stroke_alpha = 0.6,text_size = 9)+ ggtitle(paste0(xid,'_ggvenn_mrna_up.pdf'))
    #pdf(paste0(outpath,'/',xid,'_ggvenn_mrna_up.pdf'),width = 8)
    print(p)
    #dev.off()
    
    p = ggvenn::ggvenn(data = list(mrna.down.human = tmpmrna$human_downmrna,
                                   mrna.down.macaca = tmpmrna$macaca_downmrna,
                                   mrna.down.mouse = tmpmrna$mouse_downmrna),
                       fill_color = c("#E64B35FF", "#3C5488FF", "green"),
                       show_percentage = F,stroke_size = 0.5,
                       stroke_alpha = 0.6,text_size = 9) + ggtitle(paste0(xid,'_ggvenn_mrna_down.pdf'))
    #pdf(paste0(outpath,'/',xid,'_ggvenn_mrna_down.pdf'),width = 8)
    print(p)
    #dev.off()
    
    k = k+1
    # rho up
    gotype = c("C2")
    tt = 0
    for(j in 1:length(gotype)){
      
      k1 = tmpGO$mouse_upmrna[tmpGO$mouse_upmrna$type == gotype[j],]
      k2 = tmpGO$macaca_upmrna[tmpGO$macaca_upmrna$type == gotype[j],]
      k3 = tmpGO$human_upmrna[tmpGO$human_upmrna$type == gotype[j],]
      
      tt = tt+1
      vid= intersect(rownames(k1),rownames(k2))
      out.rho.up[tt,k] = cor.test(log2(k1[vid,]$ratio),log2(k2[vid,]$ratio))$estimate
      
      tt = tt+1
      vid= intersect(rownames(k2),rownames(k3))
      out.rho.up[tt,k] = cor.test(log2(k2[vid,]$ratio),log2(k3[vid,]$ratio))$estimate
      
      tt = tt+1
      vid= intersect(rownames(k1),rownames(k3))
      out.rho.up[tt,k] = cor.test(log2(k1[vid,]$ratio),log2(k3[vid,]$ratio))$estimate
    }
    
    ## go down
    # rho up
    tt = 0
    for(j in 1:length(gotype)){
      
      k1 = tmpGO$mouse_downmrna[tmpGO$mouse_downmrna$type == gotype[j],]
      k2 = tmpGO$macaca_downmrna[tmpGO$macaca_downmrna$type == gotype[j],]
      k3 = tmpGO$human_downmrna[tmpGO$human_downmrna$type == gotype[j],]
      
      tt = tt+1
      vid= intersect(rownames(k1),rownames(k2))
      out.rho.down[tt,k] = cor.test(log2(k1[vid,]$ratio),log2(k2[vid,]$ratio))$estimate
      
      tt = tt+1
      vid= intersect(rownames(k2),rownames(k3))
      out.rho.down[tt,k] = cor.test(log2(k2[vid,]$ratio),log2(k3[vid,]$ratio))$estimate
      
      tt = tt+1
      vid= intersect(rownames(k1),rownames(k3))
      out.rho.down[tt,k] = cor.test(log2(k1[vid,]$ratio),log2(k3[vid,]$ratio))$estimate
    }
   
  }
  out.rho = list(up = out.rho.up,down = out.rho.down)
  dev.off()
  return(out.rho)
}




# gene info
# mouse
mousegeneInfo = file2frame('./data/mouse_geneInfo_GRCm39v113.txt')
mousegeneInfo = mousegeneInfo[!duplicated(mousegeneInfo$Symbol) & !is.na(mousegeneInfo$Symbol),]
rownames(mousegeneInfo) = mousegeneInfo$Symbol

#human
geneInfo.v38 = file2frame('./data/geneInfo_encodev38.tab')
geneInfo.v38 = geneInfo.v38[!duplicated(geneInfo.v38$Symbol),]
rownames(geneInfo.v38) = geneInfo.v38$Symbol

#macaca
macacageneInfo = file2frame('./data/Macaca_mulatta.Mmul_10.99.ensemble_symbol_biotype.txt')
macacageneInfo = macacageneInfo[!duplicated(macacageneInfo$Symbol) & !is.na(macacageneInfo$Symbol),]
rownames(macacageneInfo) = macacageneInfo$Symbol

load('./data/rnaData.Rdata')
rnaData = updateRNA1(rnaData,headers,macacageneInfo)

mrna.tissues = rnaData$mrna.tissues
mrna.tissues.info = rnaData$mrna.tissues.info
mrna.whole = rnaData$mrna.whole
mrna.whole.info = rnaData$mrna.whole.info

#macaca
DEmrna.tissues.lm = get_tissue_DEgenes_lm(mrna.tissues,mrna.tissues.info,tissue.systems = NULL)
Metamrna = DEmrna.tissues.lm$MetaLimma

#mouse
msdata = loadRData('./data/Mus_agingNature_reads_list.Rdata')
sum(sapply(msdata$Mus_agingNature_pdata_list, nrow))
msdata$cpm = list()
msdata$cpm.clin = list()
for(i in 1:length(msdata$Mus_agingNature_reads_list)){
  counts = msdata$Mus_agingNature_reads_list[[i]]
  rownames(counts) = toupper(rownames(counts))
  clin = msdata$Mus_agingNature_pdata_list[[i]]
  rownames(clin) = clin$`Sample name`
    #consider coding genes
  tmptype = mousegeneInfo[rownames(counts),]$Type
  idx = tmptype == "protein_coding" | substr(tmptype,nchar(tmptype)-4, nchar(tmptype)) == "_gene"
  idx[is.na(idx)] = FALSE
  counts = counts[idx,]
  clin = clin[colnames(counts),]
  cpm = get_CPM_counts(counts)
  cpm.nonormal = get_CPM_nonormal(counts)
  #filtering
  idx = rowMeans(cpm.nonormal) > 1 & rowSums(counts < 5) < ncol(counts)*0.2
  cpm = cpm[idx,]
  
  clin$age = as.numeric(clin$`characteristics: age`)
  clin$type = clin$age # just for data hase the type field
  idx = clin$age >=3 & clin$`characteristics: sex` == 'f'
  cpm = cpm[,idx]
  clin = clin[idx,]
  clin$type = clin$age
  
  # del outliers
  pcamrna = prcomp(t(cpm),cor = F)
  vid = !is.outliner(pcamrna$x[,1])
  
  msdata$cpm[[i]] = cpm[,vid]
  msdata$cpm.clin[[i]] = as.data.frame(clin[vid,])
}
names(msdata$cpm) = names(msdata$Mus_agingNature_reads_list)
names(msdata$cpm.clin) = names(msdata$Mus_agingNature_reads_list)

#mean(sapply(msdata$cpm, nrow))

MetageneMouse = get_tissue_DEgenes_lm(msdata$cpm,msdata$cpm.clin,tissue.systems = NULL)$MetaLimma


# human 
#load('/home/ligh/projects/X/projects/Cell_aging/GTEX/')

load('./data/GETX.cpm_21_Female_tissues.RData')
#GTEX.cpm.v  and GTEX.clin.v



MetageneGETX = get_tissue_DEgenes_lm(GTEX.cpm.v,GTEX.clin.v,tissue.systems = NULL)$MetaLimma

#
library(msigdbr)
c2_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2")
m_t2g <- c2_gene_sets %>% dplyr::select(gs_name, gene_symbol)



outpath = './results/compare_ms_mcc_human'
slid_mrnas  = get_slid_mrnas_3_species_beta(MetageneMouse,Metamrna,MetageneGETX)
slid_GO = get_GO_slid_mrnas_1(slid_mrnas,m_t2g) # may take 10~20 mins
venfile = paste0(outpath,'GO_compare_rank_by_beta_C2enrich_venn.pdf')
slid_rho = plot_GOenrich_ratio_1(slid_mrnas,slid_GO,venfile)

pdf(paste0(outpath,'GO_compare_rank_by_beta_C2enrich_rho.pdf'),width = 4.2, height = 4.7 )

tmpdata = reshape2::melt(slid_rho$up)
tmpdata$topnumber = as.numeric(gsub("TOP","",tmpdata$Var2))

idx = substr(tmpdata$Var1,1,2) == 'C2'
ggplot(tmpdata[idx,],aes(x = topnumber, y = value, color = Var1)) + geom_line()+ geom_point()+ 
   lghplot.addtheme(size = 12,legend.position = "bottom")+ ggtitle("Up-regulated genes")+
   xlab("Number of top ranked genes") + ylab("Rho of log2ratio")

tmpdata = reshape2::melt(slid_rho$down)
tmpdata$topnumber = as.numeric(gsub("TOP","",tmpdata$Var2))

idx = substr(tmpdata$Var1,1,2) == 'C2'
ggplot(tmpdata[idx,],aes(x = topnumber, y = value, color = Var1)) + geom_line()+ geom_point()+ 
   lghplot.addtheme(size = 12,legend.position = "bottom")+ ggtitle("Down-regulated genes")+ 
   xlab("Number of top ranked genes") + ylab("Rho of log2ratio")

dev.off()


