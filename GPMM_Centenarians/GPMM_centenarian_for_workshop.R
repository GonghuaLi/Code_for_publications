options(warn = -1)    
library(ggplot2)
library(reshape2)
library(Rmisc)
library(Hmisc)
library(grid)
library(EnhancedVolcano)
library(pracma)
library(car)
library(gridExtra)
setwd("E:/github/Code_for_publications/GPMM_Centenarians")
source('./subroutines.R')

Recon3.annote = file2frame('./data/Recon3_rxns_curated.txt')
rownames(Recon3.annote) = Recon3.annote$rxns
HPM = as.matrix(file2frame('./data/human_moped_Kim_nM.txt',row.names = 1))
RNA = as.matrix(file2frame('./data/human_RNA_matrix.txt',row.names = 1))
colnames(HPM) = capitalize(colnames(HPM))
colnames(RNA) = capitalize(colnames(RNA))

# Number of EC
id_Kcat = Recon3.annote$Kcat >0;
sum(id_Kcat)
id_EC = Recon3.annote$EC_number != ''
sum(id_EC)
ratio_noKcat = sum(!id_Kcat & id_EC)/sum(id_EC)
ratio_noKcat

# plot Genes have EC but no Kcat
idcc = !id_Kcat & id_EC;
genes_noKcat = Recon3.annote$genes[idcc]
genes_noKcat = paste0(genes_noKcat,collapse = ";");
genes_noKcat = toupper(unique(strsplit(genes_noKcat,';')[[1]]));

HPM_noKcat = HPM[intersect(genes_noKcat,rownames(HPM)),]
perc_noKcat = colSums(HPM_noKcat)/colSums(HPM)

tmp = c('Ratio_no_kcat',capitalize(colnames(HPM)))
x_axis = factor(tmp,levels = tmp)

Type = c('#Reactions without Kcats', rep(c('Protein abundance without Kcats'),ncol(HPM)))
pdf(file = "./figures/for_pub/v2/Figure S1.pdf", bg = "transparent")
ggplot(,aes(x = x_axis,y = c(ratio_noKcat,perc_noKcat),color = Type,fill = Type)) + 
       geom_bar(position="dodge", stat="identity") + theme_bw()+ #ylim(c(0,1))+#scale_x_log10()+
       theme(legend.position="top") +
        xlab("Tissues")+ylab("Fraction of reaction without Kcats ") +
        lghplot.addthemeA(hjust = 1,size = 16,sizex = 12,sizey = 12,legend.position = 'top')
dev.off()

overlapgenes = intersect(rownames(RNA),rownames(HPM))
RNA.v = log2(RNA[overlapgenes,]+1)
HPM.v = log2(HPM[overlapgenes,]+1)
p =list()
for(i in 1:ncol(HPM.v)){
    idaa = RNA.v[,i] >0 & HPM.v[,i] >0
    tcor = cor.test(RNA.v[idaa,i],HPM.v[idaa,i])$estimate
    p[[i]] = ggplot(,aes(RNA.v[,i],HPM.v[,i]))+geom_point(size =1,) + theme_bw()+
        theme(plot.margin = margin(0.1,0.1,0.1,0.1,"cm"))+
        lghplot.addthemeA(size = 16,sizex = 16,sizey = 16)+
        
        annotate(geom="text", x=4, y=17, 
                label = paste('R=',signif(tcor,3)),
               color="darkblue",size = 8,face = "italic")+
        xlab('')+
        ylab('')+
         ggtitle(colnames(HPM.v)[i])
}

pdf(file = "./figures/for_pub/v2/Figure S2.pdf",)
grid.arrange(arrangeGrob(grobs = p,ncol = 4,
                         bottom=textGrob('mRNA expression(log2 FPKM)', gp=gpar(fontface="bold",  fontsize=22)),
                        left = textGrob('Protein abundance(log2 HPM)', gp=gpar(fontface="bold",  fontsize=22),rot=90)))
dev.off()


png(file = "./figures/for_pub/v2/Figure S2.png",width = 1000,height = 1000)
grid.arrange(arrangeGrob(grobs = p,ncol = 4,
                         bottom=textGrob('mRNA expression(log2 FPKM)', gp=gpar(fontface="bold",  fontsize=22)),
                        left = textGrob('Protein abundance(log2 HPM)', gp=gpar(fontface="bold",  fontsize=22),rot=90)))
dev.off()

# get ratio
HPM.overlap = HPM[overlapgenes,]
RNA.overlap = RNA[overlapgenes,]
ratio = HPM.overlap/RNA.overlap
ratio[ratio == 0] = NA
ratio[is.infinite(ratio)] = NA
ratio = apply(ratio,1,median,na.rm =T)

# remove ratio = 0 and NA
id = is.na(ratio) | ratio == 0
ratio = ratio[!id]

HPM.prediction = RNA.overlap[!id,]*ratio

pred.log2 =  log2(HPM.prediction+1)
meas.log2 = log2(HPM.overlap[!id,]+1)

p1 =list()
for(i in 1:ncol(meas.log2)){
    idaa = pred.log2[,i] > 0 & meas.log2[,i] > 0
    tcor = cor.test(pred.log2[idaa,i],meas.log2[idaa,i])$estimate
    p1[[i]] = ggplot(,aes(pred.log2[,i],meas.log2[,i]))+geom_point(size =1,) + theme_bw()+
        theme(plot.margin = margin(0.1,0.1,0.1,0.1,"cm"))+
        lghplot.addthemeA(size = 16,sizex = 16,sizey = 16)+
        
        annotate(geom="text", x=5, y=17, 
                label = paste('R=',signif(tcor,3)),
               color="darkblue",size = 8,face = "italic")+
        xlab('')+
        ylab('')+
         ggtitle(colnames(meas.log2)[i])
}

pdf(file = "./figures/for_pub/v2/Figure S3.pdf",)
grid.arrange(arrangeGrob(grobs = p1,ncol = 4,
                         bottom=textGrob('Predicted protein abundance (log2 HPM)', gp=gpar(fontface="bold",  fontsize=22)),
                        left = textGrob('Protein abundance (log2 HPM)', gp=gpar(fontface="bold",  fontsize=22),rot=90)))

dev.off()

png(file = "./figures/for_pub/v2/Figure S3.png",width = 1000,height = 1000)
grid.arrange(arrangeGrob(grobs = p1,ncol = 4,
                         bottom=textGrob('Predicted protein abundance (log2 HPM)', gp=gpar(fontface="bold",  fontsize=22)),
                        left = textGrob('Protein abundance (log2 HPM)', gp=gpar(fontface="bold",  fontsize=22),rot=90)))
dev.off()

# robust all
GPMM= as.matrix(file2frame("./benchmark/robust_all_GPMM_fluxRxnsMean.txt",row.names = 1))
flux3.annote = file2frame("./benchmark/Recon3_rxns_v1b.txt",row.names = 1)
clin = file2frame("./benchmark/NCI60_clinical59.txt")
rownames(clin) = clin$ids
benchflux = as.matrix(file2frame("./benchmark/robust_NCI60_flux_benchmark_Recon3(mmol L min)59.txt",row.names = 1))
# filter GPMM
GPMM = GPMM[rowMeans(abs(GPMM)) > 1e-6,] # remain the flux with the overall precison of >1e-6
vid = GPMM['DM_atp_c_',] > 0 # obtain valid modeling
GPMM = GPMM[,vid]
clin.GPMM = clin[vid,]
benchflux.GPMM = benchflux[,vid]


idxaa = regexpr('pert',colnames(GPMM)) < 0
GPMM.wt = GPMM[,idxaa]
corMatrix = matrix(NA,ncol(GPMM.wt),3)
rownames(corMatrix) = colnames(GPMM.wt)
perc = c('0.01','0.05','0.1')
colnames(corMatrix) = paste0('perc_',perc)
for(i in 1:length(colnames(GPMM.wt))){
    tname = colnames(GPMM.wt)[i]
    twt = GPMM.wt[,tname]
    for(j in 1:length(perc)){
        mutname = paste0(tname,'_pert',perc[j])
        if (!any(is.element(colnames(GPMM),mutname))) {next;}
        mut = GPMM[,mutname]
        corMatrix[i,j] = cor(log10(abs(twt)),log10(abs(mut)))  
    }
}

colnames(corMatrix) = c('1%','5%','10%')
xx = melt(corMatrix);

#graphics.off()
pdf("./figures/for_pub/v2/Figure 2B_robust_all.pdf")
ggplot(xx,aes(x = Var2,y = value^2))+ geom_boxplot()+lghplot.addtheme(size = 28,legend.position = 'top')+ geom_jitter(width = 0.2,size = 2)+
  xlab('NCI 60 Noised sample type')+ ylab("R-Squared")+ ylim(0.95,1)
dev.off()


#robust analysis_ 0.05noise_LC_NCI_H460
robust =as.matrix(file2frame("./benchmark/robust_GPMM_0.05noise_LC_NCI_H460_fluxRxnsMean.txt",row.names = 1))
robust = robust[rowMeans(abs(robust)) > 1e-6,] # remain the flux with the overall precison of >1e-6
vid = robust['DM_atp_c_',] > 0 # obtain valid modeling
robust = robust[,vid]

vx = intersect(rownames(GPMM),rownames(robust))
robust = cbind(GPMM[vx,'LC_NCI_H460'],robust[vx,])

#idx = abs(robust[,1]) < 10
tcor = as.vector(cor(log10(abs(robust[,1])),log10(abs(robust[,-1]))))
pdf("./figures/for_pub/v2/Figure 2C_LC_NCI_H460.pdf")

ggplot(,aes(x = 1:100,y = tcor^2)) + geom_point(color = "black")+ geom_line(color = "black")+
      lghplot.addtheme(size = 24)+  labs(title = "Robust")+ylim(0.95,1.0)+
      xlab("Noised Sample(0.05)")+ylab("R-Squared")
dev.off()

tmp = data.frame(lactate_atp = c(c(as.vector(robust['EX_lac_L(e)',])),
                               c(as.vector(robust['DM_atp_c_',]))),
                 Reaction = c(rep('Lactate secretion',101),rep('ATP production',101)),
                 stringsAsFactors = F,
                 class = c(c('WT',rep('Noised',100)),c('WT',rep('Noised',100))),
                 tx = c(1:101,1:101)
                 )
pdf("./figures/for_pub/v2/Figure 2D_robust_lactateATP.pdf")
ggplot(tmp,aes(x = tx,y = lactate_atp,color = Reaction)) + geom_point(aes(shape = class))+
      lghplot.addtheme(size = 24,legend.position = 'top')+ ylim(1.4,2.6)+
      ggsci::scale_color_aaas()+
      xlab("Noised Sample(0.05)")+ylab("mmol/min/L")
dev.off()

# read original data
GPMM= as.matrix(file2frame("./benchmark/GPMM_fluxRxnsMean.txt",row.names = 1))
GIMME = as.matrix(file2frame("./benchmark/GIMME_fluxRxnsMean.txt",row.names = 1))
FastCore = as.matrix(file2frame("./benchmark/FastCore_fluxRxnsMean.txt",row.names = 1))
rFASTCORMICS = as.matrix(file2frame("./benchmark/rFastcorm_fluxRxnsMean.txt",row.names = 1))
flux3.annote = file2frame("./benchmark/Recon3_rxns_v1b.txt",row.names = 1)
clin = file2frame("./benchmark/NCI60_clinical59.txt")
rownames(clin) = clin$ids
benchflux = as.matrix(file2frame("./benchmark/NCI60_flux_benchmark_Recon3(mmol L min)59.txt",row.names = 1))

# filter GPMM
GPMM = GPMM[rowMeans(abs(GPMM)) > 1e-6,] # remain the flux with the overall precison of >1e-6
vid = GPMM['DM_atp_c_',] > 0 # obtain valid modeling
GPMM = GPMM[,vid]
clin.GPMM = clin[vid,]
benchflux.GPMM = benchflux[,vid]

# filter GIMME
GIMME = GIMME[rowMeans(abs(GIMME)) > 1e-6,] # remain the flux with the overall precison of >1e-6
vid = GIMME['DM_atp_c_',] > 0 # obtain valid modeling
GIMME = GIMME[,vid]
clin.GIMME = clin[vid,]
benchflux.GIMME = benchflux[,vid]

# filter FastCore
FastCore = FastCore[rowMeans(abs(FastCore)) > 1e-6,] # remain the flux with the overall precison of >1e-6
vid = FastCore['DM_atp_c_',] > 0 # obtain valid modeling
FastCore = FastCore[,vid]
clin.FastCore = clin[vid,]
benchflux.FastCore = benchflux[,vid]


# filter rFASTCORMICS
rFASTCORMICS = rFASTCORMICS[rowMeans(abs(rFASTCORMICS)) > 1e-6,] # remain the flux with the overall precison of >1e-6
vid = rFASTCORMICS['DM_atp_c_',] > 0 # obtain valid modeling
rFASTCORMICS = rFASTCORMICS[,vid]
clin.rFASTCORMICS = clin[vid,]
benchflux.rFASTCORMICS = benchflux[,vid]



vrxns = intersect(rownames(benchflux.GPMM),rownames(GPMM))
prediction = GPMM[vrxns,]
experiment = benchflux.GPMM[vrxns,]

aa = apply(experiment,1,median) > 1e-3
prediction = prediction[aa,]
experiment = experiment[aa,]


for(i in 1:nrow(experiment)){
    dx =is.outliner(experiment[i,],coef = 1.5)
    experiment[i,dx] = NA
}
p.GPMM = melt(prediction)$value
e.GPMM = melt(experiment)$value

# Figure S4A
ids = e.GPMM > 0   & abs(p.GPMM) > 0
ids[is.na(ids)] = FALSE
Nrxns = sum(ids)
Nrxns
GPMMcor = cor.test(log10(abs(p.GPMM[ids])),log10(abs(e.GPMM[ids])),method = 'pearson')
GPMMcor
text = paste0('R-Squred = ',signif(GPMMcor$estimate^2,2),'\n','p = ',signif(GPMMcor$p.value,2))
text
pdf("./figures/for_pub/v2/Figure S4A.pdf")
ggplot(,aes(log10(p.GPMM[ids]),log10(e.GPMM[ids]))) + geom_point(color = "black")+ #geom_smooth(method = 'lm')+
      annotate(geom="text", x=-3, y=0.5, parse = TRUE,
               label = expression(atop(paste(italic(R)^2,' = 0.72'),
                                      paste(italic(P), '= 2.3e-106'))),
               color="darkblue",size = 12)+
      geom_abline(intercept=0,slope=1,size = 1.5,color = "darkblue")+#theme_bw()+
      annotate(geom="text", x= 0.85, y=1.35, label="1:1",color="darkblue",size = 14)+ #geom_smooth()+#ylim(-3,1.4)+
      lghplot.addtheme(size = 24)+  labs(title = "GPMM")+xlim(-5,1.35)+ ylim(-5,1.35)+
      xlab("Predicted flux (log10 mmol/min/L)")+ylab("Experimental flux (log10 mmol/min/L)")
dev.off()


vrxns = intersect(rownames(benchflux.GIMME),rownames(GIMME))
prediction = GIMME[vrxns,]
experiment = benchflux.GIMME[vrxns,]
aa = apply(experiment,1,median) > 1e-3
prediction = prediction[aa,]
experiment = experiment[aa,]

for(i in 1:nrow(experiment)){
    dx =is.outliner(experiment[i,],coef = 1.5)
    experiment[i,dx] = NA
}

p.GIMME = melt(prediction)$value
e.GIMME = melt(experiment)$value

# Figure S4B
ids = e.GIMME > 0  & abs(p.GIMME) > 0
ids[is.na(ids)] = FALSE
Nrxns = sum(ids)
Nrxns
GIMMEcor = cor.test(log10(p.GIMME[ids]),log10(e.GIMME[ids]),method = 'pearson')
GIMMEcor 
text = paste('r = ',signif(GIMMEcor$estimate^2,3),'\n','p = ',signif(GIMMEcor$p.value,3),sep = '')
text
pdf("./figures/for_pub/v2/Figure S4B.pdf")
ggplot(,aes(log10(p.GIMME[ids]),log10(e.GIMME[ids]))) + geom_point(color = "black")+ #geom_smooth(method = 'loess')+
      annotate(geom="text", x=-2, y=0.5, 
               label = expression(atop(paste(italic(R)^2,' = 0.011'),
                                      paste(italic(P), '= 0.047'))),,
               color="darkblue",size = 12,face = "italic")+
      geom_abline(intercept=0,slope=1,size = 1,color = "darkblue")+#theme_bw()+
      annotate(geom="text", x=0.6, y=1, label="1:1",color="darkblue",size = 14)+#ylim(-7,1.1)
      lghplot.addtheme(size = 24)+ labs(title = "GIMME")+xlim(-5,1)+ ylim(-5,1)+
      xlab("Predicted flux (log10 mmol/min/L)")+ylab("Experimental flux (log10 mmol/min/L)")
dev.off()


vrxns = intersect(rownames(benchflux.FastCore),rownames(FastCore))
prediction = FastCore[vrxns,]
experiment = benchflux.FastCore[vrxns,]

aa = apply(experiment,1,median) > 1e-3
prediction = prediction[aa,]
experiment = experiment[aa,]


for(i in 1:nrow(experiment)){
    dx =is.outliner(experiment[i,],coef = 1.5)
    experiment[i,dx] = NA
}

p.FastCore = melt(prediction)$value
e.FastCore = melt(experiment)$value

# Figure S4C
ids = e.FastCore > 0   & abs(p.FastCore) > 0
ids[is.na(ids)] = FALSE
Nrxns = sum(ids)
Nrxns
FastCorecor = cor.test(log10(p.FastCore[ids]),log10(e.FastCore[ids]),method = 'pearson')
FastCorecor
FastCorecor$p.value
text = paste('r = ',signif(FastCorecor$estimate^2,3),'\n','p = ',signif(FastCorecor$p.value,3),sep = '')
pdf("./figures/for_pub/v2/Figure S4C.pdf")
text
ggplot(,aes(log10(p.FastCore[ids]),log10(e.FastCore[ids]))) + geom_point(color = "black")+ #geom_smooth(method = 'loess')+
      annotate(geom="text", x=-2, y=1.36, 
               label = expression(atop(paste(italic(R)^2,' = 0.31'),
                                      paste(italic(P), '= 3.09e-24'))),
               color="darkblue",size = 12)+
      geom_abline(intercept=0,slope=1,size = 1.2,color = "darkblue")+#theme_bw()+
      annotate(geom="text", x=0.95, y=1.5, label="1:1",color="black",size = 14)+
      lghplot.addtheme(size = 24)+ labs(title = "FastCore")+xlim(-5,1.7)+ylim(-5,1.7)+
      xlab("Predicted flux (log10 mmol/min/L)")+ylab("Experimental flux (log10 mmol/min/L)")
dev.off()

vrxns = intersect(rownames(benchflux.rFASTCORMICS),rownames(rFASTCORMICS))
prediction = rFASTCORMICS[vrxns,]
experiment = benchflux.rFASTCORMICS[vrxns,]

aa = apply(experiment,1,median) > 1e-3
prediction = prediction[aa,]
experiment = experiment[aa,]


for(i in 1:nrow(experiment)){
    dx =is.outliner(experiment[i,],coef = 1.5)
    experiment[i,dx] = NA
}

p.rFASTCORMICS = melt(prediction)$value
e.rFASTCORMICS = melt(experiment)$value

# Figure 1D
ids = e.rFASTCORMICS > 0   & abs(p.rFASTCORMICS) > 0
ids[is.na(ids)] = FALSE
Nrxns = sum(ids)
Nrxns
rFASTCORMICScor = cor.test(log10(p.rFASTCORMICS[ids]),log10(e.rFASTCORMICS[ids]),method = 'pearson')
rFASTCORMICScor
rFASTCORMICScor$p.value
text = paste('r = ',signif(rFASTCORMICScor$estimate^2,3),'\n','p = ',signif(rFASTCORMICScor$p.value,3),sep = '')
pdf("./figures/for_pub/v2/Figure S4D.pdf")
text
ggplot(,aes(log10(p.rFASTCORMICS[ids]),log10(e.rFASTCORMICS[ids]))) + geom_point(color = "black")+ #geom_smooth(method = 'loess')+
      annotate(geom="text", x=-2, y=1.36, 
               label = expression(atop(paste(italic(R)^2,' = 0.488'),
                                      paste(italic(P), '= 2.01e-54'))),
               color="darkblue",size = 12)+
      geom_abline(intercept=0,slope=1,size = 1.2,color = "darkblue")+#theme_bw()+
      annotate(geom="text", x=0.95, y=1.5, label="1:1",color="black",size = 14)+
      lghplot.addtheme(size = 24)+ labs(title = "rFASTCORMICS")+xlim(-5,1.7)+ylim(-5,1.7)+
      xlab("Predicted flux (log10 mmol/min/L)")+ylab("Experimental flux (log10 mmol/min/L)")
dev.off()

ida = 'EX_lac_L(e)'
GPMMcor.lac = cor.test(GPMM[ida,],benchflux.GPMM[ida,],method = 'pearson')
#text = paste('r = ',signif(GPMMcor.lac$estimate,3),'\n','p = ',signif(GPMMcor.lac$p.value,3),sep = '')

#GPMMcor.lac

text = paste('r2 = ',signif(GPMMcor.lac$estimate^2,2),'\n','p = ',signif(GPMMcor.lac$p.value,2),sep = '')
text
pdf("./figures/for_pub/v2/Figure 2E.pdf")
ggplot(,aes(GPMM[ida,],benchflux.GPMM[ida,])) + geom_point(size = 4,color = "black")+ 

      annotate(geom="text", x=2.5, y=6, #fontface =2,#label=text,
               label = expression(atop(paste(italic(R)^2,' = 0.86'),
                                      paste(italic(P), '= 2.2e-24'))),
               
               color="darkblue",size = 12)+
      geom_abline(intercept=0,slope=1,size = 2,color = "darkblue")+#theme_bw()+
      annotate(geom="text", x=6.5, y=7, label="1:1",color="black",size = 14)+
      lghplot.addtheme(size = 26)+ ylim(0,7)+labs(title = "GPMM")+
      xlab("Predicted flux (mmol/min/L)")+ylab("Experimental flux (mmol/min/L)")
dev.off()

ida = 'EX_lac_L(e)'
GIMMEcor.lac = cor.test(GIMME[ida,],benchflux.GIMME[ida,],method = 'pearson')
text = paste('r = ',signif(GIMMEcor.lac$estimate,3),'\n','p = ',signif(GIMMEcor.lac$p.value,3),sep = '')
pdf("./figures/for_pub/v2/Figure 2F.pdf")
ggplot(,aes(GIMME[ida,],benchflux.GIMME[ida,])) + geom_point(size=4,color = "black")+ 
      #annotate(geom="text", x=2.5, y=6, label=text,color="darkblue",size = 12)+
      geom_abline(intercept=0,slope=1,size = 2,color = "darkblue")+#theme_bw()+
      annotate(geom="text", x=5.8, y=6.4, label="1:1",color="black",size = 14)+
      lghplot.addtheme(size = 26)+ ylim(0,7)+ xlim(0,6)+labs(title = "GIMME")+
      xlab("Predicted flux (mmol/min/L)")+ylab("Experimental flux (mmol/min/L)")
dev.off()

ida = 'EX_lac_L(e)'
FastCorecor.lac = cor.test(FastCore[ida,],benchflux.FastCore[ida,],method = 'pearson')
#text = paste('r = ',signif(FastCorecor.lac$estimate,3),'\n','p = ',signif(FastCorecor.lac$p.value,3),sep = '')
text = paste('r2 = ',signif(FastCorecor.lac$estimate^2,2),'\n','p = ',signif(FastCorecor.lac$p.value,2),sep = '')
text


pdf("./figures/for_pub/v2/Figure 2G.pdf")
ggplot(,aes(FastCore[ida,],benchflux.FastCore[ida,])) + geom_point(size = 4,color = "black")+ 
      annotate(geom="text", x=2.5, y=6, #label=text,
               label = expression(atop(paste(italic(R)^2,' = 0.088'),
                                      paste(italic(P), '= 0.022'))),
               color="darkblue",size = 12)+
      geom_abline(intercept=0,slope=1,size = 2,color = "darkblue")+#theme_bw()+
      annotate(geom="text", x=6.5, y=7, label="1:1",color="black",size = 14)+
      lghplot.addtheme(size = 26)+ ylim(0,7)+labs(title = "FastCore")+
      xlab("Predicted flux (mmol/min/L)")+ylab("Experimental flux (mmol/min/L)")
dev.off()

ida = 'EX_lac_L(e)'
rFASTCORMICScor.lac = cor.test(rFASTCORMICS[ida,],benchflux.rFASTCORMICS[ida,],method = 'pearson')
text = paste('r2 = ',signif(rFASTCORMICScor.lac$estimate^2,2),'\n','p = ',signif(rFASTCORMICScor.lac$p.value,2),sep = '')
text


pdf("./figures/for_pub/v2/Figure 2H.pdf")
ggplot(,aes(rFASTCORMICS[ida,],benchflux.rFASTCORMICS[ida,])) + geom_point(size = 4,color = "black")+ 
      annotate(geom="text", x=2.5, y=6, #label=text,
               label = expression(atop(paste(italic(R)^2,' = 0.33'),
                                      paste(italic(P), '= 1.8e-6'))),
               color="darkblue",size = 12)+
      geom_abline(intercept=0,slope=1,size = 2,color = "darkblue")+#theme_bw()+
      annotate(geom="text", x=6.5, y=7, label="1:1",color="black",size = 14)+
      lghplot.addtheme(size = 26)+ ylim(0,7)+labs(title = "rFASTCORMICS")+
      xlab("Predicted flux (mmol/min/L)")+ylab("Experimental flux (mmol/min/L)")
dev.off()

ecfluxes = file2frame('./benchmark/compare_ecModel/ec_GEMs/Results/11_cellLines_NCI60/ecModels_const_0_exchangeFluxesComp.txt',row.names = 1)
ecfluxes = ecfluxes[,-1]
ia = 2*(1:(ncol(ecfluxes)/2))
fluxes.ecmodel = data.frame(Allexperiment  = 7.14*melt(ecfluxes[,ia-1])$value,
                            Allprediction = 7.14*melt(ecfluxes[,ia])$value)
sum(fluxes.ecmodel$Allexperiment > 0  & fluxes.ecmodel$Allprediction < 0)
sum(fluxes.ecmodel$Allexperiment > 0  & fluxes.ecmodel$Allprediction > 0)

ids = fluxes.ecmodel$Allexperiment > 1e-3  
Nrxns = sum(ids)
SGMMcor = cor.test(log10(abs(fluxes.ecmodel$Allprediction[ids])+1e-6),log10(abs(fluxes.ecmodel$Allexperiment[ids])+1e-6),method = 'pearson')
text = paste('r = ',signif(SGMMcor$estimate^2,3),'\n','p = ',signif(SGMMcor$p.value,3),sep = '')

pdf("./figures/for_pub/Figure S4E_ecModel_overall.pdf")
text
tmpdata = data.frame(tprediction = log10(abs(fluxes.ecmodel$Allprediction[ids])+1e-6),stringsAsFactors = F,
                     tmeasured = log10(fluxes.ecmodel$Allexperiment[ids]),
                     directionType = fluxes.ecmodel$Allprediction[ids] > 0)

ggplot(tmpdata,aes(tprediction,tmeasured)) + geom_point(size = 4)+ 
      annotate(geom="text", x=-3, y=1, 
               label = expression(atop(paste(italic(R)^2,' = 0.27'),
                                      paste(italic(P), '= 3.7e-06'))),
               color="darkblue",size = 12)+
      geom_abline(intercept=0,slope=1,size = 1.2,color = "darkblue")+#theme_bw()+
      annotate(geom="text", x= 0.85, y=1.35, label="1:1",color="darkblue",size = 14)+ #geom_smooth()+#ylim(-3,1.4)+
      lghplot.addtheme(size = 24,legend.position = 'top')+  xlim(-6.1,2) +ylim(-3.1,2)+ labs(title = "ecmodel")+
      theme(legend.text=element_text(size=20,))+
      xlab("Predicted flux (log10 abs mmol/min/L)")+ylab("Experimental flux (log10 abs mmol/min/L)")
      #geom_smooth(se = FALSE, method = "gam", formula = y~x ) 
dev.off()


ecfluxes = file2frame('./benchmark/compare_ecModel/ec_GEMs/Results/11_cellLines_NCI60/ecModels_const_0_exchangeFluxesComp.txt',row.names = 1)
ecfluxes = ecfluxes['HMR_9135',-1]

ia = 2*(1:(ncol(ecfluxes)/2))
fluxes.ecmodel = data.frame(Allexperiment  = 7.14*melt(ecfluxes[,ia-1])$value,
                            Allprediction = 7.14*melt(ecfluxes[,ia])$value)

FastCorecor.lac = cor.test(fluxes.ecmodel$Allexperiment,fluxes.ecmodel$Allprediction,method = 'pearson')
#text = paste('r = ',signif(FastCorecor.lac$estimate,3),'\n','p = ',signif(FastCorecor.lac$p.value,3),sep = '')
text = paste('r2 = ',signif(FastCorecor.lac$estimate^2,2),'\n','p = ',signif(FastCorecor.lac$p.value,2),sep = '')
text


pdf("./figures/for_pub/v2/Figure 2I.pdf")
ggplot(fluxes.ecmodel,aes(Allprediction,Allexperiment)) + geom_point(size = 4,color = "black")+ 
      geom_abline(intercept=0,slope=1,size = 2,color = "darkblue")+#theme_bw()+
      annotate(geom="text", x=95, y=100, label="1:1",color="black",size = 14)+
      lghplot.addtheme(size = 26)+ labs(title = "ecModel")+xlim(0,100)+ylim(0,100)+
      xlab("Predicted flux (mmol/min/L)")+ylab("Experimental flux (mmol/min/L)")
dev.off()

flux3 = file2frame('./metabolic_modeling_centenarians/Hainan_CMWN_fluxRxnsMean_cen_fpkm_counts_cds_cmwn_combat_20191112.txt',row.names = 1)
flux3.annote = file2frame('./metabolic_modeling_centenarians/Recon3_rxns_v1b.txt',row.names = 1)
clin = file2frame('./metabolic_modeling_centenarians/Hainan_CMWN_Clinical_full.txt',row.names = 1)
list[IA,IB] = ismember(rownames(clin),colnames(flux3))
clin = clin[IA,]
flux3 = flux3[,IB]
idxx = rowMeans(abs(flux3))> 1e-6
flux3 = flux3[idxx,]
flux3.log2= as.matrix(log2(abs(flux3)+1e-6))
list[IA,IB] = ismember(rownames(flux3.log2),rownames(flux3.annote))
flux3.annote = flux3.annote[IB,]


# two linear model
# For centenarian effect:   flux ~ cen + sex     subset: not F1
# For ageing effect:  flux~ age +sex  subset: not centenariains

nrxns = nrow(flux3.log2)
lmflux3 = data.frame(rxns = rownames(flux3.log2),stringsAsFactors = F,
                    meanRxns = rowMeans(flux3.log2),
                    beta.cen = rep(0,nrxns),p.cen = rep(1,nrxns),
                    beta.age = rep(0,nrxns),p.age = rep(1,nrxns),
                    beta.gender = rep(0,nrxns),p.gender = rep(1,nrxns))
for (i in 1:nrow(flux3.log2)){
    bx = lm(flux3.log2[i,] ~ (clin$Class=='C') +  (clin$Gender == 'Female'),subset = clin$Class != 'F1')
    #bx = lm(flux3.log2[i,] ~ (clin$Class=='C') +  (clin$Class=='F1') + (clin$Gender == 'Female'))
    lmflux3$beta.cen[i] = bx$coefficients[2]
    lmflux3$p.cen[i] = car::Anova(bx)$`Pr(>F)`[1]
    #lmflux3$beta.herit[i] = bx$coefficients[3]
    #lmflux3$p.herit[i] = Anova(bx)$`Pr(>F)`[2]
    lmflux3$beta.gender[i] = bx$coefficients[3]
    lmflux3$p.gender[i] = car::Anova(bx)$`Pr(>F)`[2]
    
    bx1 = lm(flux3.log2[i,] ~ clin$Age +  (clin$Gender == 'Female'),subset = clin$Class == 'F1SP')
    ax1 = car::Anova(bx1)
    lmflux3$beta.age[i] = bx1$coefficients[2]
    lmflux3$p.age[i] = ax1$`Pr(>F)`[1]
    lmflux3$beta.gender1[i] = bx1$coefficients[3]
    lmflux3$p.gender1[i] = ax1$`Pr(>F)`[2]
}
rownames(lmflux3) = rownames(flux3.log2)
lmflux3 = cbind(lmflux3,flux3.annote)

idaa = lmflux3$p.cen < 0.05 & lmflux3$p.age < 0.05 & lmflux3$beta.cen * lmflux3$beta.age > 0
lmflux3.filter = lmflux3
lmflux3.filter$p.cen[idaa] =1

# bootstrap without replacement for lmflux3
# two linear model random
# For centenarian effect:   flux ~ cen + sex     subset: not F1
# For ageing effect:  flux~ age +sex  subset: not centenariains
#
# -------Will spend a long time------
nboot = 1000
bootstrapLmflux3 = list()
bootstrapLmflux3[[1]] = matrix(0,nrxns,nboot)  # beta.cen
bootstrapLmflux3[[2]] = matrix(1,nrxns,nboot)  # p.cen
bootstrapLmflux3[[3]] = matrix(0,nrxns,nboot)  # beta.age
bootstrapLmflux3[[4]] = matrix(1,nrxns,nboot)  # p.age
bootstrapLmflux3[[5]] = matrix(0,nrxns,nboot)  # beta.gender
bootstrapLmflux3[[6]] = matrix(1,nrxns,nboot)  # p.gender
names(bootstrapLmflux3) = c('beta.cen','p.cen','beta.age','p.age','beta.gender','p.gender')
for (k in 1:nboot){
    #randperm
    clin$RandClass = clin$Class[randperm(nrow(clin))]
    clin$RandGender = clin$Gender[randperm(nrow(clin))]
    clin$RandAge = clin$Age[randperm(nrow(clin))]
    
    for (i in 1:nrow(flux3.log2)){
        bx = lm(flux3.log2[i,] ~ (clin$RandClass=='C') +  (clin$RandGender == 'Female'),subset = clin$RandClass != 'F1')
        bootstrapLmflux3[[1]][i,k] = bx$coefficients[2]
        bootstrapLmflux3[[2]][i,k] = car::Anova(bx)$`Pr(>F)`[1]
        
        bx1 = lm(flux3.log2[i,] ~ clin$RandAge +  (clin$RandGender == 'Female'),subset = clin$RandClass == 'F1SP')
        ax1 = car::Anova(bx1)
        bootstrapLmflux3[[3]][i,k] = bx1$coefficients[2]
        bootstrapLmflux3[[4]][i,k] = ax1$`Pr(>F)`[1]
        
        bootstrapLmflux3[[5]][i,k]  = bx$coefficients[3]
        bootstrapLmflux3[[6]][i,k]  = car::Anova(bx)$`Pr(>F)`[2]
     }
}


save(list = c('flux3.log2','lmflux3','lmflux3.filter','flux3.annote','clin','bootstrapLmflux3'),file = './metabolic_modeling_centenarians/Flux3_lm_bootstrp_v20211112.Rdata')


# we suggest save bootstrap
load('./metabolic_modeling_centenarians/Flux3_lm_bootstrp_v20211112.Rdata')


# remove age effect
sum(lmflux3$beta.cen >0 &  lmflux3$p.cen < 0.05 & lmflux3$p.age < 0.05 & lmflux3$beta.age >0)
lmflux3$p.cen[lmflux3$beta.cen >0 &  lmflux3$p.cen < 0.05 & lmflux3$p.age < 0.05 & lmflux3$beta.age >0 ] = 1
sum(lmflux3$beta.cen <0 &  lmflux3$p.cen < 0.05 & lmflux3$p.age < 0.05 & lmflux3$beta.age <0)
lmflux3$p.cen[lmflux3$beta.cen <0 &  lmflux3$p.cen < 0.05 & lmflux3$p.age < 0.05 & lmflux3$beta.age <0 ] = 1

# Volcono plot of centenarians effect
# set the base colour as 'black' 
idaa = lmflux3.filter$subSystemes == 'Exchange/demand reaction' & rowMeans(flux3) < 0 & substr(rownames(lmflux3.filter),1,2) == 'EX' & rowMeans(flux3.log2) > log2(1e-5)

res2 = lmflux3.filter[idaa,1:13]
keyvals <- rep('gray50', nrow(res2))
names(keyvals) <- rep('NS', nrow(res2))   

keyvals[which(res2$beta.cen > 0 & res2$p.cen < 0.05)] <- "Brown"    
names(keyvals)[which(res2$beta.cen > 0 & res2$p.cen < 0.05)] <- 'High'    

keyvals[which(res2$beta.cen < -0 & res2$p.cen < 0.05)] <- "darkblue"    
names(keyvals)[which(res2$beta.cen < -0 & res2$p.cen < 0.05)] <- 'Low'   
pdf('./figures/for_pub/v2/Figure 3B.pdf')
EnhancedVolcano(res2, lab = rownames(res2),   x = 'beta.cen', pCutoff = 0.05,cutoffLineCol = 'white',  colAlpha = 0.9,
                legendPosition = 'botton',
                title = 'Uptake  fluxes',border = 'full',titleLabSize = 28,
                FCcutoff = 0,cutoffLineWidth = 0,axisLabSize = 28,
                gridlines.minor = F,gridlines.major = F,
                subtitle = NULL,transcriptPointSize = -5*log10(res2$p.cen),transcriptLabSize = 8,boxedlabels = F,
                y = 'p.cen',    xlim = c(-0.6, 0.6),ylim= c(0,3), xlab = "Cen effect(beta lm) ", colCustom = keyvals)
dev.off()

# calc WDA and pvalue
tmplmflux = lmflux3.filter
tmplmflux$log2FC = lmflux3.filter$beta.cen
tmplmflux$Pvalue = lmflux3.filter$p.cen
WDA.cen.bootstrap = Flux.subsystem_WDA_bootstrap(tmplmflux,bootstrapLmflux3[[1]], bootstrapLmflux3[[2]],
                                       flux3.annote,cutoff.nrxns = 3,cutoff.pval =0.05,cutoff.fc = 0)


pdf('./figures/for_pub/v2/Figure 3D.pdf',width = 11,height =8)
idaa = substr(WDA.cen.bootstrap$subsystem,1,9) == 'Transport'
print(plot_DAscore(WDA.cen.bootstrap$DAscore[idaa],WDA.cen.bootstrap$subsystem[idaa],enlarge = 1.8,gtitle = 'Transport fluxes',
                   #lengend.position = "top", 
                   WDA.cen.bootstrap$fdr[idaa],gsize =28,gsizex = 20,gsizey = 24))
dev.off()

pdf('./figures/for_pub/v2/Figure 3C.pdf',width = 15,height =18)
idaa = substr(WDA.cen.bootstrap$subsystem,1,9) != 'Transport' & WDA.cen.bootstrap$subsystem != 'Exchange/demand reaction' & WDA.cen.bootstrap$fdr < 0.05
print(plot_DAscore(WDA.cen.bootstrap$DAscore[idaa],WDA.cen.bootstrap$subsystem[idaa],enlarge = 3,gtitle = 'Enzymatic reactions',
                   WDA.cen.bootstrap$fdr[idaa],gsize =30,gsizex = 26,gsizey = 28))
dev.off()

# Volcono plot of centenarians effect
# set the base colour as 'black' 
idaa = lmflux3.filter$subSystemes == 'Exchange/demand reaction' & rowMeans(flux3) > 0 & substr(rownames(lmflux3.filter),1,2) == 'EX' & rowMeans(flux3.log2) > log2(1e-5)

res2 = lmflux3.filter[idaa,1:13]
keyvals <- rep('gray50', nrow(res2))
names(keyvals) <- rep('NS', nrow(res2))   

keyvals[which(res2$beta.cen > 0 & res2$p.cen < 0.05)] <- "Brown"    
names(keyvals)[which(res2$beta.cen > 0 & res2$p.cen < 0.05)] <- 'High'    

keyvals[which(res2$beta.cen < -0 & res2$p.cen < 0.05)] <- "darkblue"    
names(keyvals)[which(res2$beta.cen < -0 & res2$p.cen < 0.05)] <- 'Low'   
pdf('./figures/for_pub/v2/Figure 3E.pdf')
EnhancedVolcano(res2, lab = rownames(res2),   x = 'beta.cen', pCutoff = 0.05,cutoffLineCol = 'white',  colAlpha = 0.9,
                legendPosition = 'botton',
                title = 'Secretory fluxes',border = 'full',titleLabSize = 28,
                FCcutoff = 0,cutoffLineWidth = 0,axisLabSize = 28,
                gridlines.minor = F,gridlines.major = F,
                subtitle = NULL,transcriptPointSize = -3.5*log10(res2$p.cen),transcriptLabSize = 3.5,boxedlabels = F,
                y = 'p.cen',    xlim = c(-1, 1),ylim= c(0,4), xlab = "Cen effect(beta lm) ", colCustom = keyvals)
dev.off()

#text = paste('r = ',signif(tcor$estimate,2),'\n','p = ',signif(tcor$p.value,2),sep = '')
pdf("./figures/for_pub/v2/Figure_S5A.pdf",width = 5,height = 6)
tcor = cor.test(lmflux3$beta.age, lmflux3$beta.cen)
text = paste('r = ',signif(tcor$estimate,2),'\n','p = ',signif(tcor$p.value,2),sep = '')
cor.test(lmflux3$beta.age, lmflux3$beta.cen)
ggplot(,aes(lmflux3$beta.age, lmflux3$beta.cen))+lghplot.addthemeA()+xlab('Beta age')+geom_point(alpha = 0.2)+
annotate(geom="text", x=0, y=1.5, label=text,color="darkblue",size = 8)+xlim(-0.1,0.1) + ylim(-2,2)+
geom_smooth(method = 'lm')+ylab('Beta Cen')+ ggtitle('Cen effect Vs Age effect')
dev.off()

# Overlap
require(VennDiagram)
Aging_up = rownames(lmflux3)[lmflux3$beta.age > 0 & lmflux3$p.age < 0.05]
Aging_down = rownames(lmflux3)[lmflux3$beta.age < 0 & lmflux3$p.age < 0.05]
CEN_up = rownames(lmflux3)[lmflux3$beta.cen > 0 & lmflux3$p.cen < 0.05]
CEN_down = rownames(lmflux3)[lmflux3$beta.cen < 0 & lmflux3$p.cen < 0.05]
venn.diagram(list(Aging_up=Aging_up,Aging_down = Aging_down,CEN_up = CEN_up,CEN_down=CEN_down), 
             #fill=c("red","blue",),
             fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
             alpha=c(0.5,0.5,0.5,0.5), cex=2, cat.fontface=4, cat.cex = 1.3,margin = 0.05,
             filename="./figures/for_pub/v2/Figure_S5B.tiff")

# Oxidative phosphorylation complexes
lmflux3$ID = rownames(lmflux3)
ids = c('NADH2_u10mi','FADH2ETC','CYOR_u10mi','CYOOm3i','ATPS4mi')
tnames =c('Complex I','Complex II','Complex III','Complex IV','Complex V')
bx.cen = lmflux3[ids,]
tcolor = rep('p<0.05',nrow(bx.cen))
tcolor[bx.cen$p.cen > 0.05] = 'p>0.05'
bx.cen$class = tcolor
bx.cen$ID = factor(bx.cen$ID,levels = bx.cen$ID)
bx.cen$names = tnames
ida = sort.int(bx.cen$p.cen,decreasing = T,index.return = T)$ix
bx.cen = bx.cen[ida,]
bx.cen$names = factor(bx.cen$names,levels = bx.cen$names)
beta_direct = rep('UP',nrow(bx.cen))
beta_direct[bx.cen$beta.cen > 0 ] = 'UP'
beta_direct[bx.cen$beta.cen <= 0 ] = 'Down'
bx.cen$effect = as.character(beta_direct)

pdf('./figures/for_pub/v2/Figure S6A.pdf',width = 8,height = 6)
print(plot_DAscore(bx.cen$beta.cen,bx.cen$names,bx.cen$p.cen,xlabel = 'Beta.cen',gtype = 'pvalue',
                   gsize = 22,gsizex = 22,gsizey = 22,
                   gtitle = 'OP complexes',lengend.position = 'right'))
dev.off()

# ATP production
pdf("./figures/for_pub/v2/Figure S6B.pdf",width = 5,height = 6)
text = paste('p = ',signif(lmflux3.filter['DM_atp_c_',]$p.cen,2),'\n','beta = ',signif(lmflux3.filter['DM_atp_c_',]$beta.cen,2),sep = '')
id = clin$Class != 'F1'
p= ggplot(,aes(x = clin$Class[id],y = unlist(flux3.log2['DM_atp_c_',id]), fill = clin$Class[id], color = clin$Class[id])) +
lghplot.addthemeA()+ xlab('') + ylab('Log2(mmol/min/L)') + ggtitle('ATP production') + 
annotate(geom="text", x='F1SP', y=3.8, label=text,color="darkblue",size = 8)+ylim(0,4)
print(lghplot.boxplot(p))
dev.off()

DEflux.f1vsf1sp = DEGenes.simplified(flux3.log2,catagory = clin$Class == 'F1',subset = clin$Class == 'F1' | clin$Class == 'F1SP')

dim(lmflux3)
dim(DEflux.f1vsf1sp)

pdf("./figures/for_pub/v2/Figure S9A.pdf",width = 6,height = 6)
tcor = cor.test(lmflux3$beta.cen,DEflux.f1vsf1sp$log2FC)
text = paste('r = ',signif(tcor$estimate,2),'\n','p = ',signif(tcor$p.value,2),sep = '')
cor.test(lmflux3$beta.cen,DEflux.f1vsf1sp$log2FC)
ggplot(,aes(lmflux3$beta.cen,DEflux.f1vsf1sp$log2FC))+lghplot.addthemeA()+geom_point(alpha = 0.2)+
annotate(geom="text", x=-0.5, y=0.8, label=text,color="darkblue",size = 8)+xlim(-1,1) + ylim(-1,1)+
geom_smooth(method = 'lm')+xlab('Beta Cen')+ ylab('F1 vs F1SP log2FC')+ggtitle('Beta Cen Vs F1s log2FC')
dev.off()

idx1 = lmflux3$beta.cen > 0 &  lmflux3$p.cen < 0.05 & lmflux3$subSystemes == 'Oxidative phosphorylation'
idx2 = lmflux3$beta.cen > 0 &  lmflux3$p.cen < 0.05 & lmflux3$subSystemes == 'Fatty acid beta oxidation'
idx3 = lmflux3$beta.cen > 0 &  lmflux3$p.cen < 0.05 & lmflux3$subSystemes == 'Transport, peroxisomal'
idx4 = lmflux3$beta.cen > 0 &  lmflux3$p.cen < 0.05 & lmflux3$subSystemes == 'Citric acid cycle'
idx5 = rownames(lmflux3) == 'DM_atp_c_'
#idx2 = lmflux3$beta.cen > 0 &  lmflux3$p.cen < 0.05 & lmflux3$subSystemes == 'Transport, peroxisomal'
idxx = idx1 | idx2 | idx3 | idx4 | idx5
sum(idxx)
sum(DEflux.f1vsf1sp$log2FC[idxx] >0)
#plot(DEflux.f1vsf1sp$log2FC[idxx],lmflux3$beta.cen[idxx])

#tmpdata = cbind(lmflux3[idxx,],DEflux.f1vsf1sp[idxx,])
res2 = cbind(lmflux3[idxx,],DEflux.f1vsf1sp[idxx,])
res2 = res2[,c('beta.cen','log2FC','Pvalue')]
#res2$beta.cen = 2^res2$beta.cen
keyvals <- rep('gray50', nrow(res2))
names(keyvals) <- rep('NS', nrow(res2))   

keyvals[which(res2$log2FC > 0 )] <- "Brown"    
names(keyvals)[which(res2$log2FC > 0 )] <- 'High'    

keyvals[which(res2$log2FC < -0 )] <- "darkblue"    
names(keyvals)[which(res2$log2FC < -0 )] <- 'Low'
res2$keyvals = keyvals
pdf("./figures/for_pub/v2/Figure S9B.pdf",width = 6,height = 6)
ggplot(res2,aes(x = log2FC, y = beta.cen)) + 
    geom_point(color = keyvals,aes(size = beta.cen*3))+
    geom_text_repel(label = rownames(res2)) + theme_bw()+lghplot.addtheme()+ 
    xlab('F1 vs F1SP log2FC')+ ylab('Beta Cen') +ggtitle('Beta Cen Vs F1s log2F:FAO')
dev.off()




#writetxt(lmflux3[,c('ID','beta.cen')],'./metabolic_modeling_centenarians/Hainan_CMWN_Flux3_lm_result_20191112_for_Escher.txt',sep = ',',row.names= F, col.names = F)

# Figure 2G
#  draw Figures in EScher:   https://escher.github.io
#  step 1: open https://escher.github.io and select map "Glycolysis TCA PPP " with the Model "none"
#  step 2: uplaod Recon3 model of "Recon3v1_for_Escher.json" in "metabolic_modeling_centenarians" subdir
#  step 3:  upload Reaction fluxes data of "Hainan_CMWN_Flux3_lm_result_20191112_for_Escher.txt" in "metabolic_modeling_centenarians" subdir
#  step 4: draw Figure 2G  in escher.

expr = as.matrix(file2frame('./metabolic_modeling_centenarians/GPMM_for_CEN/data/cen_fpkm_counts_cds_cmwn_combat_20191112.txt',row.names = 1))
expr = log2(expr[,colnames(flux3)]+1)

idx = flux3.annote$subSystemes == 'Fatty acid beta oxidation'
faogenes = toupper(unique(unlist(strsplit(unique(flux3.annote$genes[idx]),';'))))
expr.fao = expr[intersect(faogenes,rownames(expr)),]



nrxns = nrow(expr.fao)
lmfaogene = data.frame(expr = rownames(expr.fao),stringsAsFactors = F,
                    meanExpr = rowMeans(expr.fao),
                    beta.cen = rep(0,nrxns),p.cen = rep(1,nrxns),
                    beta.age = rep(0,nrxns),p.age = rep(1,nrxns),
                    beta.gender = rep(0,nrxns),p.gender = rep(1,nrxns))
for (i in 1:nrow(expr.fao)){
    bx = lm(expr.fao[i,] ~ (clin$Class=='C') +  (clin$Gender == 'Female'),subset = clin$Class != 'F1')
    #bx = lm(flux3.log2[i,] ~ (clin$Class=='C') +  (clin$Class=='F1') + (clin$Gender == 'Female'))
    lmfaogene$beta.cen[i] = bx$coefficients[2]
    lmfaogene$p.cen[i] = car::Anova(bx)$`Pr(>F)`[1]
    #lmflux3$beta.herit[i] = bx$coefficients[3]
    #lmflux3$p.herit[i] = Anova(bx)$`Pr(>F)`[2]
    
    bx1 = lm(expr.fao[i,] ~ clin$Age +  (clin$Gender == 'Female'),subset = clin$Class != 'C')
    ax1 = car::Anova(bx1)
    lmfaogene$beta.age[i] = bx1$coefficients[2]
    lmfaogene$p.age[i] = ax1$`Pr(>F)`[1]
    lmfaogene$beta.gender[i] = bx1$coefficients[3]
    lmfaogene$p.gender[i] = ax1$`Pr(>F)`[2]
}
rownames(lmfaogene) = rownames(expr.fao)
#lmfaogene = cbind(lmfaogene,flux3.annote)


#cen
res2 = lmfaogene
keyvals <- rep('gray50', nrow(res2))
names(keyvals) <- rep('NS', nrow(res2))   

keyvals[which(res2$beta.cen > 0 & res2$p.cen < 0.05)] <- "Brown"    
names(keyvals)[which(res2$beta.cen > 0 & res2$p.cen < 0.05)] <- 'High'    

keyvals[which(res2$beta.cen < -0 & res2$p.cen < 0.05)] <- "darkblue"    
names(keyvals)[which(res2$beta.cen < -0 & res2$p.cen < 0.05)] <- 'Low' 
tlab = rownames(res2)
pdf("./figures/for_pub/v2/Figure S8A.pdf")
EnhancedVolcano(res2, lab = tlab,   x = 'beta.cen', pCutoff = 0.4,cutoffLineCol = 'white', colAlpha = 0.9,#title = 'Beta FAO related gene expression change',
                legendPosition = 'botton',
                border = 'full',titleLabSize = 28,
                FCcutoff = 0,cutoffLineWidth = 0,axisLabSize = 28,
                gridlines.minor = F,gridlines.major = F,
                ylim = c(0,1.05*max(-log10(res2$p.cen), na.rm=TRUE)),xlim= c(-0.1,0.16),
                subtitle = NULL,transcriptPointSize = -3*log10(res2$p.cen),transcriptLabSize = 6,boxedlabels = F,
                y = 'p.cen',  xlab = "Cen effect(beta lm)", colCustom = keyvals)
dev.off()

#age
res2 = lmfaogene
keyvals <- rep('gray50', nrow(res2))
names(keyvals) <- rep('NS', nrow(res2))   

keyvals[which(res2$beta.age > 0 & res2$p.age < 0.05)] <- "Brown"    
names(keyvals)[which(res2$beta.age > 0 & res2$p.age < 0.05)] <- 'High'    

keyvals[which(res2$beta.age < -0 & res2$p.age < 0.05)] <- "darkblue"    
names(keyvals)[which(res2$beta.age < -0 & res2$p.age < 0.05)] <- 'Low' 
tlab = rownames(res2)
pdf("./figures/for_pub/v2/Figure S8B.pdf")
EnhancedVolcano(res2, lab = tlab,   x = 'beta.age', pCutoff = 0.4,cutoffLineCol = 'white', colAlpha = 0.9,#title = 'Beta FAO related gene expression change',
                legendPosition = 'botton',
                border = 'full',titleLabSize = 28,
                FCcutoff = 0,cutoffLineWidth = 0,axisLabSize = 28,
                gridlines.minor = F,gridlines.major = F,
                ylim = c(0,1.05*max(-log10(res2$p.age), na.rm=TRUE)),xlim= c(-0.011,0.007),
                subtitle = NULL,transcriptPointSize = -3*log10(res2$p.age),transcriptLabSize = 6,boxedlabels = F,
                y = 'p.age',  xlab = "Age effect(beta lm)", colCustom = keyvals)
dev.off()

met = file2frame('./metabolic_modeling_centenarians/Hainan_CMWN_metabolism_addhmdb-filter0.8.txt')
clin.met = file2frame('./metabolic_modeling_centenarians/Hainan_CMWN_Clinical_full.txt',row.names = 1)
rownames(met) = paste0('Met_',1:nrow(met))
met.header = met[,c(1:17)]
met.expr = fillgaps_rowMin(as.matrix(met[,-c(1:17)]))
met.expr = log2(met.expr)
list[IA,IB] = ismember(colnames(met.expr),rownames(clin.met))
met.expr = met.expr[,IA]
clin.met = clin.met[IB,]

dim(met.expr)

nrxns = nrow(met.expr)
lmMet = data.frame(meanRxns = rowMeans(met.expr),
                    beta.cen = rep(0,nrxns),p.cen = rep(1,nrxns),
                    beta.age = rep(0,nrxns),p.age = rep(1,nrxns),
                    beta.gender = rep(0,nrxns),p.gender = rep(1,nrxns))
for (i in 1:nrow(met.expr)){
   bx = lm(met.expr[i,] ~ (clin.met$Class=='C') +  (clin.met$Gender == 'Female'),subset = clin.met$Class != 'F1')
    #bx = lm(met.expr[i,] ~ (clin.met$Class=='C') +  (clin.met$Class=='F1') + (clin.met$Gender == 'Female'))
    lmMet$beta.cen[i] = bx$coefficients[2]
    lmMet$p.cen[i] = car::Anova(bx)$`Pr(>F)`[1] 
    lmMet$beta.gender[i] = bx$coefficients[3]
    lmMet$p.gender[i] = car::Anova(bx)$`Pr(>F)`[2] 
    
    bx1 = lm(met.expr[i,] ~ clin.met$Age +  (clin.met$Gender == 'Female'),subset = clin.met$Class == 'F1SP')
    ax1 = car::Anova(bx1)
    lmMet$beta.age[i] = bx1$coefficients[2]
    lmMet$p.age[i] = ax1$`Pr(>F)`[1]

}
rownames(lmMet) = rownames(met.expr)
lmMet = cbind(met.header,lmMet)
lmMet$ID = rownames(lmMet)

idaa = lmMet$p.cen < 0.05 & lmMet$p.age < 0.05 & lmMet$beta.cen * lmMet$beta.age > 0
lmMet.filter = lmMet
lmMet.filter$p.cen[idaa] =1

# bootstrap without replacement for lmMet
# two linear model random
# For centenarian effect:   flux ~ cen + sex     subset: not F1
# For ageing effect:  flux~ age +sex  subset: not centenariains
nmets = nrow(lmMet)
nboot = 1000
bootstraplmMet = list()
bootstraplmMet[[1]] = matrix(0,nmets,nboot)  # beta.cen
bootstraplmMet[[2]] = matrix(1,nmets,nboot)  # p.cen
bootstraplmMet[[3]] = matrix(0,nmets,nboot)  # beta.age
bootstraplmMet[[4]] = matrix(1,nmets,nboot)  # p.age
bootstraplmMet[[5]] = matrix(0,nmets,nboot)  # beta.gender
bootstraplmMet[[6]] = matrix(1,nmets,nboot)  # p.gender
names(bootstraplmMet) = c('beta.cen','p.cen','beta.age','p.age','beta.gender','p.gender')
for (k in 1:nboot){
    #randperm
    clin.met$RandClass = clin.met$Class[randperm(nrow(clin.met))]
    clin.met$RandGender = clin.met$Gender[randperm(nrow(clin.met))]
    clin.met$RandAge = clin.met$Age[randperm(nrow(clin.met))]
    
    for (i in 1:nrow(met.expr)){
        bx = lm(met.expr[i,] ~ (clin.met$RandClass=='C') +  (clin.met$RandGender == 'Female'),subset = clin.met$RandClass != 'F1')
        bootstraplmMet[[1]][i,k] = bx$coefficients[2]
        bootstraplmMet[[2]][i,k] = car::Anova(bx)$`Pr(>F)`[1]
        
        bx1 = lm(met.expr[i,] ~ clin.met$RandAge +  (clin.met$RandGender == 'Female'),subset = clin.met$RandClass != 'C')
        ax1 = car::Anova(bx1)
        bootstraplmMet[[3]][i,k] = bx1$coefficients[2]
        bootstraplmMet[[4]][i,k] = ax1$`Pr(>F)`[1]
        
        bootstraplmMet[[5]][i,k]  = bx$coefficients[3]
        bootstraplmMet[[6]][i,k]  = car::Anova(bx)$`Pr(>F)`[2]
     }
}


save(list = c('lmMet','lmMet.filter','met.expr','clin.met','met.header','bootstraplmMet'),file = './metabolic_modeling_centenarians/Met_lm_bootstrp__v20211112.Rdata')


#save(list = c('lmMet','met.expr','clin.met','met.header','bootstraplmMet'),file = './metabolic_modeling_centenarians/Met_lm_bootstrp.Rdata')
#load('./metabolic_modeling_centenarians/Met_lm_bootstrp.Rdata')
#
load('./metabolic_modeling_centenarians/Met_lm_bootstrp__v20211112.Rdata')
#

# cen
res2 = lmMet
keyvals <- rep('gray50', nrow(res2))
names(keyvals) <- rep('NS', nrow(res2))   

keyvals[which(res2$beta.cen > 0 & res2$p.cen < 0.05)] <- "Brown"    
names(keyvals)[which(res2$beta.cen > 0 & res2$p.cen < 0.05)] <- 'High'    

keyvals[which(res2$beta.cen < -0 & res2$p.cen < 0.05)] <- "darkblue"    
names(keyvals)[which(res2$beta.cen < -0 & res2$p.cen < 0.05)] <- 'Low'   
pdf('./figures/for_pub/v2/Figure 4A.pdf')
EnhancedVolcano(res2, lab = res2$met_name,   x = 'beta.cen', pCutoff = 0.05,cutoffLineCol = 'white', 
                legendPosition = 'botton',
                title = 'Cen plasma metabolism',border = 'full',titleLabSize = 28,
                FCcutoff = 0,cutoffLineWidth = 0,axisLabSize = 28,
                gridlines.minor = F,gridlines.major = F,
                subtitle = NULL,transcriptPointSize = -1.2*log10(res2$p.cen),transcriptLabSize = 4,boxedlabels = F,
                y = 'p.cen',    xlim = c(-2, 4),ylim= c(0,6), xlab = "Cen effect(beta lm) ", colCustom = keyvals)
dev.off()

FALclass = c('Long-chain fatty acids','Glycerophosphocholines',
            'Phosphatidic acid', 'Glycerophosphoethanolamines', 'Glycerophosphoglycerols',
            'Glycerophosphoinositols','Sphingomyelin','Fatty amides')
ida = ismember(lmMet.filter$metClass,FALclass)[[1]]

sum(lmMet.filter$p.cen < 0.05 & ida & lmMet.filter$beta.cen < 0)
tdata = data.frame(regulateType = c('Up regulated','Up regulated','Down regulated','Down regulated'),stringsAsFactors = F,
                   metaboType = c('FAL','non_FAL','FAL','non_FAL'),
                  number = c(sum(lmMet.filter$p.cen < 0.05 & ida & lmMet.filter$beta.cen > 0),
                            sum(lmMet.filter$p.cen < 0.05 & !ida & lmMet.filter$beta.cen > 0),
                            sum(lmMet.filter$p.cen < 0.05 & ida & lmMet.filter$beta.cen < 0),
                            sum(lmMet.filter$p.cen < 0.05 & !ida & lmMet.filter$beta.cen < 0))
                   )
tdata
pdf('./figures/for_pub/v2/Figure 4B.pdf')
ggplot(data=tdata, aes(x=regulateType, y=number, fill=metaboType))+geom_col(position = "fill",)+
      ylab('Relative ratio')+ xlab('')+
      annotate(geom="text", x='Down regulated', y=0.55, label='67',color="darkblue",size = 15)+
      annotate(geom="text", x='Down regulated', y=0.10, label='14',color="darkblue",size = 15)+
      annotate(geom="text", x='Up regulated', y=0.85, label='14',color="darkblue",size = 15)+
      annotate(geom="text", x='Up regulated', y=0.35, label='39',color="darkblue",size = 15)+
      lghplot.addthemeA(legend.position = 'top',size = 22,sizex = 20, sizey = 20)
dev.off()

pdf('./figures/for_pub/v2/Figure 4C.pdf',width = 13,height = 10)
lmMet$log2FC = lmMet$beta.cen
lmMet$Pvalue = lmMet$p.cen
lmMet$subSystemes = lmMet$metClass
ids = c(1:6,13,15)
DAscoreMet = Flux.subsystem_WDA_bootstrap(lmMet,bootstraplmMet[[1]],bootstraplmMet[[2]],lmMet,cutoff.nrxns = 3)
DAscoreMet.FAL = DAscoreMet[ids,]
plot_DAscore(DAscoreMet.FAL$DAscore,subsystem = DAscoreMet.FAL$subsystem,fdr = DAscoreMet.FAL$fdr,enlarge = 2.5,gtitle = 'Metabolite class',
            gsize = 28,gsizex = 22,gsizey = 24)
dev.off()

# lipid data
lipidtype = c('PA','PC','PE','SM','PI')
mtype = substr(lmMet.filter$met_name,1,2)
ida =  is.element(mtype,lipidtype) & lmMet.filter$p.cen < 0.05
DElipid = lmMet.filter[ida,]
met.expr.DElipid = met.expr[DElipid$ID,]
controlMet = rowMeans(met.expr.DElipid[,clin.met$Class == 'F1SP' & clin.met$Gender == 'Female'])
idCen = clin.met$Class == 'C' & clin.met$Gender == 'Female'
met.expr.fc = met.expr.DElipid[,idCen] -
             matrix(rep(controlMet,times = sum(idCen)),nrow = nrow(met.expr.DElipid),ncol=sum(idCen))
           
id_CF1SP = (clin.met$Class == 'C' | clin.met$Class == 'F1SP') & clin.met$Gender == 'Female'
met.expr.log2fc = met.expr.DElipid[,]


dim(met.expr.log2fc)

tmpdata = data.frame(log2fc = as.vector(unlist(met.expr.fc)),
                     met_name = rep(DElipid$met_name,times = ncol(met.expr.fc)),stringsAsFactors = F
                     )
tgc <- summarySE(tmpdata, measurevar="log2fc", groupvars="met_name")
tgc$tcolor = rep("red",nrow(tgc))
tgc$tcolor[tgc$log2fc < 0] = "darkblue"
tgc$ymin = tgc$log2fc-tgc$se
tgc$ymax = tgc$log2fc+tgc$se
tgc$ymin[tgc$log2fc >0] = tgc$log2fc[tgc$log2fc >0]
tgc$ymax[tgc$log2fc <0] = tgc$log2fc[tgc$log2fc <0]

#sorting
lipidtype = c('PA','PC','PE','PI','SM')
tgc1 = data.frame()
for (i in 1:length(lipidtype)){
    idx = which(substr(tgc$met_name,1,2) == lipidtype[i])
    tmpx = tgc[idx,]
    tmpx = tmpx[sort.int(tmpx$log2fc,decreasing = T,index.return = T)$ix,]
    tgc1 = rbind(tgc1,tmpx)
}
tgc1$met_name = factor(tgc1$met_name,levels = tgc1$met_name)
# plot
pdf("./figures/for_pub/v2/Figure 4D.pdf",width = 4,height = 8)

ggplot(tgc1, aes(x=met_name, y=log2fc)) + 
    geom_bar(aes(fill = tcolor),position=position_dodge(0.8),color = 'gray10',stat="identity",alpha = 0.7) + 
    geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.8,size = 0.8,alpha = 0.7,  
                  position=position_dodge(.9))+theme_bw() +coord_flip()+                  
lghplot.addtheme(hjust = F,size =12)+ xlab("") +  labs(y=expression(Log[2]("Fold change"))) +scale_fill_manual(values=c('darkblue','Brown'))
dev.off()


pdf('./figures/for_pub/v2/Figure 4E.pdf')
id = lmMet$met_name == 'trans-Vaccenic acid' & lmMet$score > 0.8
ida = clin.met$Class != 'gea'
p = ggplot(,aes(x = clin.met$Class[ida],y = met.expr[id,ida],color = clin.met$Class[ida],fill = clin.met$Class[ida]))+
lghplot.addthemeA(sizex = 24,sizey = 22,size = 22) + labs(y=expression(Log[2]("Peak Area")))  + xlab('')+labs(title = "Trans-Vaccenic")
print(lghplot.boxplot(p))
dev.off()
pdf('./figures/for_pub/v2/Figure 4F.pdf')
id = lmMet$met_name == 'Palmitic acid' & lmMet$score > 0.8
ida = clin.met$Class != 'gea'
p = ggplot(,aes(x = clin.met$Class[ida],y = met.expr[id,ida],color = clin.met$Class[ida],fill = clin.met$Class[ida]))+
lghplot.addthemeA(sizex = 24,sizey = 22,size = 22) + labs(y=expression(Log[2]("Peak Area")))  + xlab('')+labs(title = "Palmitic acid")
print(lghplot.boxplot(p))
dev.off()

#pvalue
id = lmMet$met_name == 'trans-Vaccenic acid' & lmMet$score > 0.8
ida = clin.met$Class != 'F1'
t.test(met.expr[id,ida]~clin.met$Class[ida])


id = lmMet$met_name == 'trans-Vaccenic acid' & lmMet$score > 0.8
ida = clin.met$Class != 'C'
t.test(met.expr[id,ida]~clin.met$Class[ida])


id = lmMet$met_name == 'Palmitic acid' & lmMet$score > 0.8
ida = clin.met$Class != 'F1'
t.test(met.expr[id,ida]~clin.met$Class[ida])

id = lmMet$met_name == 'Palmitic acid' & lmMet$score > 0.8
ida = clin.met$Class != 'C'
t.test(met.expr[id,ida]~clin.met$Class[ida])

DEmet.f1vsf1sp = DEGenes.simplified(met.expr,catagory = clin.met$Class == 'F1',subset = clin.met$Class == 'F1SP' | clin.met$Class == 'F1')
DEmet.f1vsf1sp = cbind(lmMet[,1:17],DEmet.f1vsf1sp)

FALclass = c('Long-chain fatty acids','Glycerophosphocholines',
            'Phosphatidic acid', 'Glycerophosphoethanolamines', 'Glycerophosphoglycerols',
            'Glycerophosphoinositols','Sphingomyelin','Fatty amides')
ida = ismember(DEmet.f1vsf1sp$metClass,FALclass)[[1]]

sum(DEmet.f1vsf1sp$Pvalue < 0.05 & ida & DEmet.f1vsf1sp$log2FC < 0)
tdata = data.frame(regulateType = c('Up regulated','Up regulated','Down regulated','Down regulated'),stringsAsFactors = F,
                   metaboType = c('FAL','non_FAL','FAL','non_FAL'),
                  number = c(sum(DEmet.f1vsf1sp$Pvalue < 0.05 & ida & DEmet.f1vsf1sp$log2FC > 0),
                            sum(DEmet.f1vsf1sp$Pvalue < 0.05 & !ida & DEmet.f1vsf1sp$log2FC > 0),
                            sum(DEmet.f1vsf1sp$Pvalue < 0.05 & ida & DEmet.f1vsf1sp$log2FC < 0),
                            sum(DEmet.f1vsf1sp$Pvalue < 0.05 & !ida & DEmet.f1vsf1sp$log2FC < 0))
                   )
tdata
pdf('./figures/for_pub/v2/Figure S9C.pdf',width = 6,height = 6)
ggplot(data=tdata, aes(x=regulateType, y=number, fill=metaboType))+geom_col(position = "fill",)+
      ylab('Relative ratio')+ xlab('')+
      annotate(geom="text", x='Down regulated', y=0.55, label='12',color="darkblue",size = 15)+
      annotate(geom="text", x='Down regulated', y=0.10, label='8',color="darkblue",size = 15)+
      annotate(geom="text", x='Up regulated', y=0.85, label='7',color="darkblue",size = 15)+
      annotate(geom="text", x='Up regulated', y=0.35, label='24',color="darkblue",size = 15)+
      lghplot.addthemeA(legend.position = 'top',size = 22,sizex = 20, sizey = 20)
dev.off()



#### clin.met 
bx = lm(clin.met$TC ~ (clin.met$Class=='C') +  (clin.met$Gender == 'Female'),subset = clin.met$Class != 'F1')
tpval = car::Anova(bx)$`Pr(>F)`[1]
tbeta = bx$coefficients[2]

pdf("./figures/for_pub/v2/Figure S7.pdf")
text = paste('p = ',signif(tpval,2),'\n','beta = ',signif(tbeta,2),sep = '')
id = clin.met$Class != 'F1'
p= ggplot(,aes(x = clin.met$Class[id],y = clin.met$TC[id], fill = clin.met$Class[id], color = clin.met$Class[id])) +
lghplot.addtheme()+ xlab('') + ylab('TC (umol/L)') + ggtitle('Serum total cholesterol') + 
annotate(geom="text", x='C', y=9, label=text,color="darkblue",size = 12)
print(lghplot.boxplot(p))
dev.off()

graphics.off()
