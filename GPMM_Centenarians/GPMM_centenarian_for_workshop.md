
## 1 set pwd and subroutines


```R
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
```
    

## 2 Design GPMM

### 2.1 Read data


```R
Recon3.annote = file2frame('./data/Recon3_rxns_curated.txt')
rownames(Recon3.annote) = Recon3.annote$rxns
HPM = as.matrix(file2frame('./data/human_moped_Kim_nM.txt',row.names = 1))
RNA = as.matrix(file2frame('./data/human_RNA_matrix.txt',row.names = 1))
colnames(HPM) = capitalize(colnames(HPM))
colnames(RNA) = capitalize(colnames(RNA))
```

### 2.2 Kcat related Figure S2


```R
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
pdf(file = "./figures/for_pub/Figure S2.pdf", bg = "transparent")
ggplot(,aes(x = x_axis,y = c(ratio_noKcat,perc_noKcat),color = Type,fill = Type)) + 
       #geom_smooth(level=0.9) + theme_bw()+ #scale_x_log10()+
       geom_bar(position="dodge", stat="identity") + theme_bw()+ #ylim(c(0,1))+#scale_x_log10()+
       theme(legend.position="top") +
       #guides(colour=FALSE,fill=FALSE)+
        xlab("Tissues")+ylab("Fraction of reaction without Kcats ") +
        lghplot.addthemeA(hjust = 1,size = 16,sizex = 12,sizey = 12,legend.position = 'top')
dev.off()
```


2602



4517



0.424839495240204



<strong>png:</strong> 2


### 2.3 mRNA vs Protein abundance: Figure S3


```R
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
```


```R
pdf(file = "./figures/for_pub/Figure S3.pdf",)
grid.arrange(arrangeGrob(grobs = p,ncol = 4,
                         bottom=textGrob('mRNA expression(log2 FPKM)', gp=gpar(fontface="bold",  fontsize=22)),
                         #bottom = 'mRNA expression(log2 FPKM)',
                        left = textGrob('Protein abundance(log2 HPM)', gp=gpar(fontface="bold",  fontsize=22),rot=90)))
dev.off()

```


<strong>png:</strong> 2



```R
png(file = "./figures/for_pub/Figure S3.png",width = 1000,height = 1000)
grid.arrange(arrangeGrob(grobs = p,ncol = 4,
                         bottom=textGrob('mRNA expression(log2 FPKM)', gp=gpar(fontface="bold",  fontsize=22)),
                         #bottom = 'mRNA expression(log2 FPKM)',
                        left = textGrob('Protein abundance(log2 HPM)', gp=gpar(fontface="bold",  fontsize=22),rot=90)))
dev.off()
```


<strong>png:</strong> 2


### 2.4 predicted Protein vs Measured Protein: Figure S4


```R
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

```


```R
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
```


```R
pdf(file = "./figures/for_pub/Figure S4.pdf",)
grid.arrange(arrangeGrob(grobs = p1,ncol = 4,
                         bottom=textGrob('Predicted protein abundance (log2 HPM)', gp=gpar(fontface="bold",  fontsize=22)),
                         #bottom = 'mRNA expression(log2 FPKM)',
                        left = textGrob('Protein abundance (log2 HPM)', gp=gpar(fontface="bold",  fontsize=22),rot=90)))

dev.off()
```


<strong>png:</strong> 2



```R
png(file = "./figures/for_pub/Figure S4.png",width = 1000,height = 1000)
grid.arrange(arrangeGrob(grobs = p1,ncol = 4,
                         bottom=textGrob('Predicted protein abundance (log2 HPM)', gp=gpar(fontface="bold",  fontsize=22)),
                         #bottom = 'mRNA expression(log2 FPKM)',
                        left = textGrob('Protein abundance (log2 HPM)', gp=gpar(fontface="bold",  fontsize=22),rot=90)))
dev.off()
```


<strong>png:</strong> 2


## 3 Benchmarking GPMM

### 3.1 Read data


```R
# read original data
GPMM= as.matrix(file2frame("./benchmark/GPMM_fluxRxnsMean.txt",row.names = 1))
GIMME = as.matrix(file2frame("./benchmark/GIMME_fluxRxnsMean.txt",row.names = 1))
FastCore = as.matrix(file2frame("./benchmark/FastCore_fluxRxnsMean.txt",row.names = 1))
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

```

### 3.2 Benchmark GPMM: Figure 1B


```R
vrxns = intersect(rownames(benchflux.GPMM),rownames(GPMM))
prediction = GPMM[vrxns,]
experiment = benchflux.GPMM[vrxns,]
p.GPMM = melt(prediction)$value
e.GPMM = melt(experiment)$value

# Figure 1B
ids = e.GPMM > 1e-3   & abs(p.GPMM) > 1e-6
Nrxns = sum(ids)
Nrxns
GPMMcor = cor.test(log10(abs(p.GPMM[ids])),log10(abs(e.GPMM[ids])),method = 'pearson')
text = paste0('R-Squred = ',signif(GPMMcor$estimate^2,2),'\n','p = ',signif(GPMMcor$p.value,2))
#text = paste0(expression(R^2),' = ',signif(GPMMcor$estimate^2,2),'\n','p = ',signif(GPMMcor$p.value,2))
#text = expression(paste0('R'^2,' = ',signif(GPMMcor$estimate^2,2),'\n','p = ',signif(GPMMcor$p.value,2))
text
#expression(italic(R)^2,'= 0.65\n',italic(P),'= 1.9e-100'),
pdf("./figures/for_pub/Figure 1B.pdf")
ggplot(,aes(log10(p.GPMM[ids]),log10(e.GPMM[ids]))) + geom_point(color = "black")+ 
      annotate(geom="text", x=-3, y=1, parse = TRUE,
               #label = "italic(R) '= 25.60'  italic(p)<0.001",
                #label  = "italic(R)^2' = 0.65'~ ~ italic(p)<0.001",
               label = expression(atop(paste(italic(R)^2,' = 0.65'),
                                      paste(italic(P), '= 1.9e-100'))),
               #label = text,
               #label=expression(paste(italic(R)^2,'= 0.65 ')),
               color="darkblue",size = 12)+
      geom_abline(intercept=0,slope=1,size = 1.2,color = "darkblue")+#theme_bw()+
      annotate(geom="text", x= 0.85, y=1.35, label="1:1",color="darkblue",size = 14)+ #geom_smooth()+#ylim(-3,1.4)+
      lghplot.addtheme(size = 24)+  labs(title = "GPMM")+
      xlab("Predicted flux (log10 mmol/min/L)")+ylab("Experimental flux (log10 mmol/min/L)")
      #geom_smooth(se = FALSE, method = "gam", formula = y~x ) 
dev.off()
```


437



'R-Squred = 0.65\np = 1.9e-100'



<strong>png:</strong> 2


### 3.3 Benchmark GIMME: Figure 1C


```R
vrxns = intersect(rownames(benchflux.GIMME),rownames(GIMME))
prediction = GIMME[vrxns,]
experiment = benchflux.GIMME[vrxns,]
p.GIMME = melt(prediction)$value
e.GIMME = melt(experiment)$value

# Figure 1C
ids = e.GIMME > 1e-3   & abs(p.GIMME) > 1e-6
Nrxns = sum(ids)
Nrxns
GIMMEcor = cor.test(log10(p.GIMME[ids]),log10(e.GIMME[ids]),method = 'pearson')
text = paste('r = ',signif(GIMMEcor$estimate^2,3),'\n','p = ',signif(GIMMEcor$p.value,3),sep = '')
text
pdf("./figures/for_pub/Figure 1C.pdf")
ggplot(,aes(log10(p.GIMME[ids]),log10(e.GIMME[ids]))) + geom_point(color = "black")+ 
      annotate(geom="text", x=-2, y=0.5, 
               label = expression(atop(paste(italic(R)^2,' = 0.0018'),
                                      paste(italic(P), '= 0.40'))),,
               color="darkblue",size = 12,face = "italic")+
      geom_abline(intercept=0,slope=1,size = 1.2,color = "darkblue")+#theme_bw()+
      annotate(geom="text", x=0.1, y=0.5, label="1:1",color="darkblue",size = 14)+#ylim(-7,1.1)
      lghplot.addtheme(size = 24)+ labs(title = "GIMME")+
      xlab("Predicted flux (log10 mmol/min/L)")+ylab("Experimental flux (log10 mmol/min/L)")
dev.off()

```


385



'r = 0.00184\np = 0.401'



<strong>png:</strong> 2


### 3.4 Benchmark FastCore: Figure 1D


```R
vrxns = intersect(rownames(benchflux.FastCore),rownames(FastCore))
prediction = FastCore[vrxns,]
experiment = benchflux.FastCore[vrxns,]
p.FastCore = melt(prediction)$value
e.FastCore = melt(experiment)$value

# Figure 1D
ids = e.FastCore > 1e-3   & abs(p.FastCore) > 1e-6
Nrxns = sum(ids)
FastCorecor = cor.test(log10(p.FastCore[ids]),log10(e.FastCore[ids]),method = 'pearson')
FastCorecor
text = paste('r = ',signif(FastCorecor$estimate^2,3),'\n','p = ',signif(FastCorecor$p.value,3),sep = '')
pdf("./figures/for_pub/Figure 1D.pdf")
text
ggplot(,aes(log10(p.FastCore[ids]),log10(e.FastCore[ids]))) + geom_point(color = "black")+ # geom_smooth()+
      annotate(geom="text", x=-3, y=1.3, 
               label = expression(atop(paste(italic(R)^2,' = 0.26'),
                                      paste(italic(P), '= 1.1e-23'))),
               color="darkblue",size = 12)+
      geom_abline(intercept=0,slope=1,size = 1.2,color = "darkblue")+#theme_bw()+
      annotate(geom="text", x=0.85, y=1.3, label="1:1",color="black",size = 14)+ylim(-3.1,1.55) +
      lghplot.addtheme(size = 24)+ labs(title = "FastCore")+
      xlab("Predicted flux (log10 mmol/min/L)")+ylab("Experimental flux (log10 mmol/min/L)")
dev.off()
```


    
    	Pearson's product-moment correlation
    
    data:  log10(p.FastCore[ids]) and log10(e.FastCore[ids])
    t = 10.847, df = 335, p-value < 2.2e-16
    alternative hypothesis: true correlation is not equal to 0
    95 percent confidence interval:
     0.4262038 0.5848072
    sample estimates:
          cor 
    0.5098248 
    



'r = 0.26\np = 1.08e-23'



<strong>png:</strong> 2


### 3.5 Prediction lactate  GPMM: Figure 1E


```R
ida = 'EX_lac_L(e)'
GPMMcor.lac = cor.test(GPMM[ida,],benchflux.GPMM[ida,],method = 'pearson')
#text = paste('r = ',signif(GPMMcor.lac$estimate,3),'\n','p = ',signif(GPMMcor.lac$p.value,3),sep = '')

#GPMMcor.lac

text = paste('r2 = ',signif(GPMMcor.lac$estimate^2,2),'\n','p = ',signif(GPMMcor.lac$p.value,2),sep = '')
text
pdf("./figures/for_pub/Figure 1E.pdf")
ggplot(,aes(GPMM[ida,],benchflux.GPMM[ida,])) + geom_point(color = "black")+ 

      annotate(geom="text", x=2.5, y=6, #label=text,
               label = expression(atop(paste(italic(R)^2,' = 0.86'),
                                      paste(italic(P), '= 2.2e-24'))),
               
               color="darkblue",size = 12)+
      geom_abline(intercept=0,slope=1,size = 1.2,color = "darkblue")+#theme_bw()+
      annotate(geom="text", x=6.5, y=7, label="1:1",color="black",size = 14)+
      lghplot.addtheme(size = 24)+ ylim(0,7)+labs(title = "GPMM")+
      xlab("Predicted flux (mmol/min/L)")+ylab("Experimental flux (mmol/min/L)")
dev.off()
```


'r2 = 0.86\np = 2.2e-24'



<strong>png:</strong> 2


### 3.6  Prediction lactate GIMME: Figure 1F


```R
ida = 'EX_lac_L(e)'
GIMMEcor.lac = cor.test(GIMME[ida,],benchflux.GIMME[ida,],method = 'pearson')
text = paste('r = ',signif(GIMMEcor.lac$estimate,3),'\n','p = ',signif(GIMMEcor.lac$p.value,3),sep = '')
pdf("./figures/for_pub/Figure 1F.pdf")
ggplot(,aes(GIMME[ida,],benchflux.GIMME[ida,])) + geom_point(color = "black")+ 
      #annotate(geom="text", x=2.5, y=6, label=text,color="darkblue",size = 12)+
      geom_abline(intercept=0,slope=1,size = 1.2,color = "darkblue")+#theme_bw()+
      annotate(geom="text", x=5.8, y=6.4, label="1:1",color="black",size = 14)+
      lghplot.addtheme(size = 24)+ ylim(0,7)+ xlim(0,6)+labs(title = "GIMME")+
      xlab("Predicted flux (mmol/min/L)")+ylab("Experimental flux (mmol/min/L)")
dev.off()
```


<strong>png:</strong> 2


### 3.7 Prediction lactate FastCore: Figure 1G


```R
ida = 'EX_lac_L(e)'
FastCorecor.lac = cor.test(FastCore[ida,],benchflux.FastCore[ida,],method = 'pearson')
#text = paste('r = ',signif(FastCorecor.lac$estimate,3),'\n','p = ',signif(FastCorecor.lac$p.value,3),sep = '')
text = paste('r2 = ',signif(FastCorecor.lac$estimate^2,2),'\n','p = ',signif(FastCorecor.lac$p.value,2),sep = '')
text


pdf("./figures/for_pub/Figure 1G.pdf")
ggplot(,aes(FastCore[ida,],benchflux.FastCore[ida,])) + geom_point(color = "black")+ 
      annotate(geom="text", x=2.5, y=6, #label=text,
               label = expression(atop(paste(italic(R)^2,' = 0.088'),
                                      paste(italic(P), '= 0.022'))),
               color="darkblue",size = 12)+
      geom_abline(intercept=0,slope=1,size = 1.2,color = "darkblue")+#theme_bw()+
      annotate(geom="text", x=6.5, y=7, label="1:1",color="black",size = 14)+
      lghplot.addtheme(size = 24)+ ylim(0,7)+labs(title = "FastCore")+
      xlab("Predicted flux (mmol/min/L)")+ylab("Experimental flux (mmol/min/L)")
dev.off()
```


'r2 = 0.088\np = 0.022'



<strong>png:</strong> 2


### 3.8 compare ecModel from science signaling 2020

#### 3.8.1 from paper


```R
ecfluxes = file2frame('./benchmark/compare_ecModel/ec_GEMs/Results/11_cellLines_NCI60/ecModels_const_0_exchangeFluxesComp.txt',row.names = 1)
ecfluxes = ecfluxes[,-1]
ecfluxes = ecfluxes[rowMeans(abs(ecfluxes)) > 1e-6/7.14,  ] # remove low precision
ia = 2*(1:(ncol(ecfluxes)/2))
fluxes.ecmodel = data.frame(Allexperiment  = 7.14*melt(ecfluxes[,ia-1])$value,
                            Allprediction = 7.14*melt(ecfluxes[,ia])$value)
sum(fluxes.ecmodel$Allexperiment > 0  & fluxes.ecmodel$Allprediction < 0)
sum(fluxes.ecmodel$Allexperiment > 0  & fluxes.ecmodel$Allprediction > 0)

ids = fluxes.ecmodel$Allexperiment > 1e-3  #& abs(fluxes.ecmodel$Allprediction) > 1e-6
Nrxns = sum(ids)
SGMMcor = cor.test(log10(abs(fluxes.ecmodel$Allprediction[ids])+1e-6),log10(abs(fluxes.ecmodel$Allexperiment[ids])+1e-6),method = 'pearson')
#SGMMcor = cor.test(log10(abs(fluxes.ecmodel$Allprediction[ids])+1e-6),log10(abs(fluxes.ecmodel$Allexperiment[ids]+1e-6)),method = 'pearson')
text = paste('r = ',signif(SGMMcor$estimate^2,3),'\n','p = ',signif(SGMMcor$p.value,3),sep = '')
#pdf("./figures/Figure2A.pdf")
text
tmpdata = data.frame(tprediction = log10(abs(fluxes.ecmodel$Allprediction[ids])+1e-6),stringsAsFactors = F,
                     tmeasured = log10(fluxes.ecmodel$Allexperiment[ids]),
                     directionType = fluxes.ecmodel$Allprediction[ids] > 0)

ggplot(tmpdata,aes(tprediction,tmeasured,color = directionType)) + geom_point(size = 4)+ 
      annotate(geom="text", x=-3, y=1, 
               label = expression(atop(paste(italic(R)^2,' = 0.27'),
                                      paste(italic(P), '= 3.7e-06'))),
               color="darkblue",size = 12)+
      geom_abline(intercept=0,slope=1,size = 1.2,color = "darkblue")+#theme_bw()+
      annotate(geom="text", x= 0.85, y=1.35, label="1:1",color="darkblue",size = 14)+ #geom_smooth()+#ylim(-3,1.4)+
      lghplot.addtheme(size = 24,legend.position = 'top')+  xlim(-6.1,2) +ylim(-3.1,2)+ #labs(title = "ecmodel")+
      theme(legend.text=element_text(size=20,))+
      xlab("Predicted flux (log10 abs mmol/min/L)")+ylab("Experimental flux (log10 abs mmol/min/L)")
      #geom_smooth(se = FALSE, method = "gam", formula = y~x ) 
#dev.off()

```

    No id variables; using all as measure variables
    No id variables; using all as measure variables
    


43



34



'r = 0.265\np = 3.71e-06'



![png](output_33_4.png)



```R
fluxes.ecmodela = file2frame('./benchmark/compare_ecModel/Results/model_NCI60_science2012_uptake_constraint.txt')
x1 = sum(fluxes.ecmodel$Allexperiment >0 & fluxes.ecmodel$Allprediction > 0)
tdata = data.frame(regulateType = c('ecModel_pred1','ecModel_pred1',
                                    'ecModel_pred2','ecModel_pred2',
                                   'GPMM','GPMM'),stringsAsFactors = F,
                   directionType = c('True','False','True','False','True','False'),
                  number = c(sum(fluxes.ecmodel$Allexperiment >0 & fluxes.ecmodel$Allprediction > 0),
                            sum(fluxes.ecmodel$Allexperiment >0 & fluxes.ecmodel$Allprediction <= 0),
                            sum(fluxes.ecmodela$Allexperiment >0 & fluxes.ecmodela$Allprediction > 0),
                            sum(fluxes.ecmodela$Allexperiment >0 & fluxes.ecmodela$Allprediction <= 0),
                            sum(e.GPMM > 0 & p.GPMM > 0),
                            sum(e.GPMM > 0 & p.GPMM <= 0))
                   )
#pdf('./figures/ecModel_secretion_comparison_ratio.pdf')
ggplot(data=tdata, aes(x=regulateType, y=number, fill=directionType))+geom_col(position = "fill",)+
      ylab('Relative ratio')+ xlab('')+# scale_fill_manual(values=c('#A3A3A3','#FF83FA'))+
      #annotate(geom="text", x='Down regulated', y=0.55, label='60',color="darkblue",size = 15)+
      #annotate(geom="text", x='Down regulated', y=0.10, label='13',color="darkblue",size = 15)+
      #annotate(geom="text", x='Up regulated', y=0.85, label='16',color="darkblue",size = 15)+
      #annotate(geom="text", x='Up regulated', y=0.35, label='39',color="darkblue",size = 15)+
      lghplot.addtheme(legend.position = 'top')+theme(legend.text=element_text(size=20,))
#dev.off()

```


![png](output_34_0.png)


## 4 Metabolic modeling of centenarians using GPMM

### 4.1 Read data


```R
flux3 = file2frame('./metabolic_modeling_centenarians/Hainan_CMWN_fluxRxnsMean_cen_fpkm_counts_cds_cmwn_combat_20191112.txt',row.names = 1)
flux3.annote = file2frame('./metabolic_modeling_centenarians/Recon3_rxns_v1b.txt',row.names = 1)
clin = file2frame('./metabolic_modeling_centenarians/Hainan_CMWN_Clinical_simple.txt',row.names = 1)
list[IA,IB] = ismember(rownames(clin),colnames(flux3))
clin = clin[IA,]
flux3 = flux3[,IB]
idxx = rowMeans(abs(flux3))> 1e-6
flux3 = flux3[idxx,]
flux3.log2= as.matrix(log2(abs(flux3)+1e-6))
list[IA,IB] = ismember(rownames(flux3.log2),rownames(flux3.annote))
flux3.annote = flux3.annote[IB,]

```

### 4.2 Differential flux analysis


```R
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
    
    bx1 = lm(flux3.log2[i,] ~ clin$Age +  (clin$Gender == 'Female'),subset = clin$Class != 'C')
    ax1 = car::Anova(bx1)
    lmflux3$beta.age[i] = bx1$coefficients[2]
    lmflux3$p.age[i] = ax1$`Pr(>F)`[1]
    lmflux3$beta.gender[i] = bx1$coefficients[3]
    lmflux3$p.gender[i] = ax1$`Pr(>F)`[2]
}
rownames(lmflux3) = rownames(flux3.log2)
lmflux3 = cbind(lmflux3,flux3.annote)

```


```R
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
        bx1 = lm(flux3.log2[i,] ~ clin$RandAge +  (clin$RandGender == 'Female'),subset = clin$RandClass != 'C')
        ax1 = car::Anova(bx1)
        bootstrapLmflux3[[3]][i,k] = bx1$coefficients[2]
        bootstrapLmflux3[[4]][i,k] = ax1$`Pr(>F)`[1]
        bootstrapLmflux3[[5]][i,k]  = bx1$coefficients[3]
        bootstrapLmflux3[[6]][i,k]  = ax1$`Pr(>F)`[2]
     }
}

```


```R
# we suggest save bootstrap
#save(list = c('flux3.log2','lmflux3','flux3.annote','clin','bootstrapLmflux3'),file = './metabolic_modeling_centenarians/Flux3_lm_bootstrp.Rdata')
load('./metabolic_modeling_centenarians/Flux3_lm_bootstrp.Rdata')
rownames(lmflux3) = rownames(flux3.log2)
lmflux3 = cbind(lmflux3,flux3.annote)
```


```R
# remove age effect
sum(lmflux3$beta.cen >0 &  lmflux3$p.cen < 0.05 & lmflux3$p.age < 0.05 & lmflux3$beta.age >0)
lmflux3$p.cen[lmflux3$beta.cen >0 &  lmflux3$p.cen < 0.05 & lmflux3$p.age < 0.05 & lmflux3$beta.age >0 ] = 1
sum(lmflux3$beta.cen <0 &  lmflux3$p.cen < 0.05 & lmflux3$p.age < 0.05 & lmflux3$beta.age <0)
lmflux3$p.cen[lmflux3$beta.cen <0 &  lmflux3$p.cen < 0.05 & lmflux3$p.age < 0.05 & lmflux3$beta.age <0 ] = 1
```


6



0


#### 4.2.1 Uptake change : Figure 2B


```R
# centerians 
lmflux3$ID = rownames(lmflux3)
# uptake  nutritions seen in simulation uptake file
ids = c("EX_ocdcea(e)","EX_ocdca(e)","EX_hdca(e)",# fatty acids
        'EX_glc(e)',# glucose
       "EX_arg_L(e)","EX_gln_L(e)","EX_his_L(e)","EX_ile_L(e)","EX_leu_L(e)","EX_lys_L(e)","EX_met_L(e)","EX_phe_L(e)","EX_thr_L(e)","EX_trp_L(e)","EX_val_L(e)",# amino acid
       "EX_pi(e)","EX_na1(e)","EX_fe2(e)","EX_k(e)","EX_h2o(e)","EX_Tyr_ggn(e)","EX_o2(e)")
tnames = c('Octadecenoate', 'Octadecanoate', 'Hexadecanoate',
       'D-Glucose', 'L-Arginine','L-Glutamine', 'L-Histidine',
       'L-Isoleucine', 'L-Leucine', 'L-Lysine' ,'L-Methionine',
        'L-Phenylalanine', 'L-Threonine', 'L-Tryptophan', 'L-Valine', 'Phosphate',
        'Sodium', 'Iron(Fe2+)', 'Kalium', 'Water','Primer_for_Glycogen', 'Oxygen')
nutritiontype = c(rep('Fatty acid',3),'Glucose',rep('Amino acid',11),rep('Cofactors and Irons',7))
bx.cen = lmflux3[ids,]
tcolor = rep('p<0.05',nrow(bx.cen))
tcolor[bx.cen$p.cen > 0.05] = 'p>0.05'
bx.cen$class = tcolor
bx.cen$ID = factor(bx.cen$ID,levels = bx.cen$ID)
bx.cen$names = tnames
bx.cen$nutritiontype = nutritiontype
ida = sort.int(bx.cen$p.cen,decreasing = T,index.return = T)$ix
bx.cen = bx.cen[ida,]
bx.cen = rbind(bx.cen[bx.cen$nutritiontype == 'Cofactors and Irons',],bx.cen[bx.cen$nutritiontype == 'Amino acid',],
           bx.cen[bx.cen$nutritiontype == 'Glucose',],bx.cen[bx.cen$nutritiontype == 'Fatty acid',])
bx.cen$names = factor(bx.cen$names,levels = bx.cen$names)
beta_direct = rep('UP',nrow(bx.cen))
beta_direct[bx.cen$beta.cen > 0 ] = 'UP'
beta_direct[bx.cen$beta.cen <= 0 ] = 'Down'
bx.cen$effect = as.character(beta_direct)

#pdf('./figures/Figure_uptake_lmbased_cen_DAplot.pdf')
print(plot_DAscore(bx.cen$beta.cen,bx.cen$names,bx.cen$p.cen,xlabel = 'Beta.cen',gtype = 'pvalue',
                   gsize = 18,gtitle = 'Uptake change in CEN',lengend.position = 'right'))
#dev.off()
```


![png](output_44_0.png)



```R
# Volcono plot of centenarians effect
# set the base colour as 'black' 
idaa = lmflux3$subSystemes == 'Exchange/demand reaction' & rowMeans(flux3) < 0 & substr(rownames(lmflux3),1,2) == 'EX' & rowMeans(flux3.log2) > log2(1e-5)

res2 = lmflux3[idaa,1:13]
keyvals <- rep('gray50', nrow(res2))
names(keyvals) <- rep('NS', nrow(res2))   

keyvals[which(res2$beta.cen > 0 & res2$p.cen < 0.05)] <- "Brown"    
names(keyvals)[which(res2$beta.cen > 0 & res2$p.cen < 0.05)] <- 'High'    

keyvals[which(res2$beta.cen < -0 & res2$p.cen < 0.05)] <- "darkblue"    
names(keyvals)[which(res2$beta.cen < -0 & res2$p.cen < 0.05)] <- 'Low'   
pdf('./figures/for_pub/Figure 2B.pdf')
EnhancedVolcano(res2, lab = rownames(res2),   x = 'beta.cen', pCutoff = 0.05,cutoffLineCol = 'white',  colAlpha = 0.9,
                legendPosition = 'botton',
                title = 'Uptake  fluxes',border = 'full',titleLabSize = 28,
                FCcutoff = 0,cutoffLineWidth = 0,axisLabSize = 28,
                gridlines.minor = F,gridlines.major = F,
                subtitle = NULL,transcriptPointSize = -5*log10(res2$p.cen),transcriptLabSize = 8,boxedlabels = F,
                y = 'p.cen',    xlim = c(-0.6, 0.6),ylim= c(0,3), xlab = "Cen effect(beta lm) ", colCustom = keyvals)
dev.off()
```


<strong>png:</strong> 2


#### 4.2.2 Figure 2C WDA transport


```R
# calc WDA and pvalue
lmflux3$log2FC = lmflux3$beta.cen
lmflux3$Pvalue = lmflux3$p.cen
WDA.cen.bootstrap = Flux.subsystem_WDA_bootstrap(lmflux3,bootstrapLmflux3[[1]], bootstrapLmflux3[[2]],
                                       flux3.annote,cutoff.nrxns = 3,cutoff.pval =0.05,cutoff.fc = 0)

```


```R
pdf('./figures/for_pub/Figure 2C.pdf',width = 11,height =8)
idaa = substr(WDA.cen.bootstrap$subsystem,1,9) == 'Transport'
print(plot_DAscore(WDA.cen.bootstrap$DAscore[idaa],WDA.cen.bootstrap$subsystem[idaa],#lengend.position = "top",
                   WDA.cen.bootstrap$fdr[idaa],gsize =28,gsizex = 20,gsizey = 24))
dev.off()
```


<strong>png:</strong> 2


#### 4.2.3 Figure  2D DA enzymatic reaction


```R
pdf('./figures/for_pub/Figure 2D.pdf',width = 12,height =10)
idaa = substr(WDA.cen.bootstrap$subsystem,1,9) != 'Transport' & WDA.cen.bootstrap$subsystem != 'Exchange/demand reaction'
print(plot_DAscore(WDA.cen.bootstrap$DAscore[idaa],WDA.cen.bootstrap$subsystem[idaa],
                   WDA.cen.bootstrap$fdr[idaa],gsize =28,gsizex = 12,gsizey = 16))
dev.off()
```


<strong>png:</strong> 2


#### 4.2.4 Figure 2E volcono plot  of secretion


```R
# Volcono plot of centenarians effect
# set the base colour as 'black' 
idaa = lmflux3$subSystemes == 'Exchange/demand reaction' & rowMeans(flux3) > 0 & substr(rownames(lmflux3),1,2) == 'EX' & rowMeans(flux3.log2) > log2(1e-5)

res2 = lmflux3[idaa,1:13]
keyvals <- rep('gray50', nrow(res2))
names(keyvals) <- rep('NS', nrow(res2))   

keyvals[which(res2$beta.cen > 0 & res2$p.cen < 0.05)] <- "Brown"    
names(keyvals)[which(res2$beta.cen > 0 & res2$p.cen < 0.05)] <- 'High'    

keyvals[which(res2$beta.cen < -0 & res2$p.cen < 0.05)] <- "darkblue"    
names(keyvals)[which(res2$beta.cen < -0 & res2$p.cen < 0.05)] <- 'Low'   
pdf('./figures/for_pub/Figure 2E.pdf')
EnhancedVolcano(res2, lab = rownames(res2),   x = 'beta.cen', pCutoff = 0.05,cutoffLineCol = 'white',  colAlpha = 0.9,
                legendPosition = 'botton',
                title = 'Secretory fluxes',border = 'full',titleLabSize = 28,
                FCcutoff = 0,cutoffLineWidth = 0,axisLabSize = 28,
                gridlines.minor = F,gridlines.major = F,
                subtitle = NULL,transcriptPointSize = -3.5*log10(res2$p.cen),transcriptLabSize = 3.5,boxedlabels = F,
                y = 'p.cen',    xlim = c(-1, 1),ylim= c(0,4), xlab = "Cen effect(beta lm) ", colCustom = keyvals)
dev.off()
```


<strong>png:</strong> 2


#### 4.2.5 Figure 2F HMR3406 and HMR3407


```R
#id = clin$Gender == 'Female' & clin$Class != 'F1'
pdf("./figures/for_pub/Figure 2F_HMR_3406.pdf",width = 4,height = 6)
text = paste('p = ',signif(lmflux3['HMR_3406',]$p.cen,2),'\n','beta = ',signif(lmflux3['HMR_3406',]$beta.cen,2),sep = '')
id = clin$Class != 'F1'
p= ggplot(,aes(x = clin$Class[id],y = unlist(flux3.log2['HMR_3406',id]), fill = clin$Class[id], color = clin$Class[id])) +
lghplot.addtheme()+ xlab('') + ylab('Log2(mmol/min/L)') + ggtitle('HMR_3406') + 
annotate(geom="text", x='F1SP', y=-5.5, label=text,color="darkblue",size = 7)
print(lghplot.boxplot(p))
dev.off()
```


<strong>png:</strong> 2



```R
#id = clin$Gender == 'Female' & clin$Class != 'F1'
pdf("./figures/for_pub/Figure 2F_HMR_3407.pdf",width = 4,height = 6)
text = paste('p = ',signif(lmflux3['HMR_3407',]$p.cen,2),'\n','beta = ',signif(lmflux3['HMR_3407',]$beta.cen,2),sep = '')
id = clin$Class != 'F1'
p= ggplot(,aes(x = clin$Class[id],y = unlist(flux3.log2['HMR_3407',id]), fill = clin$Class[id], color = clin$Class[id])) +
lghplot.addtheme()+ xlab('') + ylab('') + ggtitle('HMR_3407')+ ylim(-6,-4)+
annotate(geom="text", x='C', y=-4.2, label=text,color="darkblue",size = 7)
print(lghplot.boxplot(p))
dev.off()
```


<strong>png:</strong> 2


#### 4.2.6 OP  complexes Figure S5


```R
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

pdf('./figures/for_pub/Figure S5A.pdf',width = 8,height = 6)
print(plot_DAscore(bx.cen$beta.cen,bx.cen$names,bx.cen$p.cen,xlabel = 'Beta.cen',gtype = 'pvalue',
                   gsize = 22,gsizex = 22,gsizey = 22,
                   gtitle = 'OP complexes',lengend.position = 'right'))
dev.off()
```


<strong>png:</strong> 2



```R
# ATP production
pdf("./figures/for_pub/Figure S5B.pdf",width = 5,height = 6)
text = paste('p = ',signif(lmflux3['DM_atp_c_',]$p.cen,2),'\n','beta = ',signif(lmflux3['DM_atp_c_',]$beta.cen,2),sep = '')
id = clin$Class != 'F1'
p= ggplot(,aes(x = clin$Class[id],y = unlist(flux3.log2['DM_atp_c_',id]), fill = clin$Class[id], color = clin$Class[id])) +
lghplot.addthemeA()+ xlab('') + ylab('Log2(mmol/min/L)') + ggtitle('ATP production') + 
annotate(geom="text", x='F1SP', y=3.8, label=text,color="darkblue",size = 8)+ylim(0,4)
print(lghplot.boxplot(p))
dev.off()
```


<strong>png:</strong> 2


### 4.3 Draw Figure 2G glycolysis-TCA-FAO


```R
#writetxt(lmflux3[,c('ID','beta.cen')],'./metabolic_modeling_centenarians/Hainan_CMWN_Flux3_lm_result_20191112_for_Escher.txt',sep = ',',row.names= F, col.names = F)
```


```R
# Figure 2G
#  draw Figures in EScher:   https://escher.github.io
#  step 1: open https://escher.github.io and select map "Glycolysis TCA PPP " with the Model "none"
#  step 2: uplaod Recon3 model of "Recon3v1_for_Escher.json" in "metabolic_modeling_centenarians" subdir
#  step 3:  upload Reaction fluxes data of "Hainan_CMWN_Flux3_lm_result_20191112_for_Escher.txt" in "metabolic_modeling_centenarians" subdir
#  step 4: draw Figure 2G  in escher.
```

## 5 Centenarians' plasma metabolic analysis

### 5.1 Read data


```R
met = file2frame('./metabolic_modeling_centenarians/Hainan_CMWN_metabolism_addhmdb-filter0.8.txt')
clin.met = file2frame('./metabolic_modeling_centenarians/Hainan_CMWN_Clinical_simple.txt',row.names = 1)
rownames(met) = paste0('Met_',1:nrow(met))
met.header = met[,c(1:17)]
met.expr = fillgaps_rowMin(as.matrix(met[,-c(1:17)]))
met.expr = log2(met.expr)
list[IA,IB] = ismember(colnames(met.expr),rownames(clin.met))
met.expr = met.expr[,IA]
clin.met = clin.met[IB,]
```

### 5.2 Differentail DEmet


```R
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
    bx1 = lm(met.expr[i,] ~ clin.met$Age +  (clin.met$Gender == 'Female'),subset = clin.met$Class != 'C')
    ax1 = car::Anova(bx1)
    lmMet$beta.age[i] = bx1$coefficients[2]
    lmMet$p.age[i] = ax1$`Pr(>F)`[1]
    lmMet$beta.gender[i] = bx1$coefficients[3]
    lmMet$p.gender[i] = ax1$`Pr(>F)`[2]
}
rownames(lmMet) = rownames(met.expr)
lmMet = cbind(met.header,lmMet)
lmMet$ID = rownames(lmMet)
```


```R
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
        bootstraplmMet[[5]][i,k]  = bx1$coefficients[3]
        bootstraplmMet[[6]][i,k]  = ax1$`Pr(>F)`[2]
     }
}

```


```R
#save(list = c('lmMet','met.expr','clin.met','met.header','bootstraplmMet'),file = './metabolic_modeling_centenarians/Met_lm_bootstrp.Rdata')
#load('./metabolic_modeling_centenarians/Met_lm_bootstrp.Rdata')
#
```

### 5.3 volcono metabolites: Figure 3A


```R
# cen
res2 = lmMet
keyvals <- rep('gray50', nrow(res2))
names(keyvals) <- rep('NS', nrow(res2))   

keyvals[which(res2$beta.cen > 0 & res2$p.cen < 0.05)] <- "Brown"    
names(keyvals)[which(res2$beta.cen > 0 & res2$p.cen < 0.05)] <- 'High'    

keyvals[which(res2$beta.cen < -0 & res2$p.cen < 0.05)] <- "darkblue"    
names(keyvals)[which(res2$beta.cen < -0 & res2$p.cen < 0.05)] <- 'Low'   
pdf('./figures/for_pub/Figure 3A.pdf')
EnhancedVolcano(res2, lab = res2$met_name,   x = 'beta.cen', pCutoff = 0.05,cutoffLineCol = 'white', 
                legendPosition = 'botton',
                title = 'Cen plasma metabolism',border = 'full',titleLabSize = 28,
                FCcutoff = 0,cutoffLineWidth = 0,axisLabSize = 28,
                gridlines.minor = F,gridlines.major = F,
                subtitle = NULL,transcriptPointSize = -1.2*log10(res2$p.cen),transcriptLabSize = 4,boxedlabels = F,
                y = 'p.cen',    xlim = c(-2, 4),ylim= c(0,6), xlab = "Cen effect(beta lm) ", colCustom = keyvals)
dev.off()
```


<strong>png:</strong> 2


### 5.4 FAL ratio Figure 3B


```R
FALclass = c('Long-chain fatty acids','Glycerophosphocholines',
            'Phosphatidic acid', 'Glycerophosphoethanolamines', 'Glycerophosphoglycerols',
            'Glycerophosphoinositols','Sphingomyelin','Fatty amides')
ida = ismember(lmMet$metClass,FALclass)[[1]]

sum(lmMet$p.cen < 0.05 & ida & lmMet$beta.cen < 0)
tdata = data.frame(regulateType = c('Up regulated','Up regulated','Down regulated','Down regulated'),stringsAsFactors = F,
                   metaboType = c('FAL','non_FAL','FAL','non_FAL'),
                  number = c(sum(lmMet$p.cen < 0.05 & ida & lmMet$beta.cen > 0),
                            sum(lmMet$p.cen < 0.05 & !ida & lmMet$beta.cen > 0),
                            sum(lmMet$p.cen < 0.05 & ida & lmMet$beta.cen < 0),
                            sum(lmMet$p.cen < 0.05 & !ida & lmMet$beta.cen < 0))
                   )
tdata
pdf('./figures/for_pub/Figure 3B.pdf')
ggplot(data=tdata, aes(x=regulateType, y=number, fill=metaboType))+geom_col(position = "fill",)+
      ylab('Relative ratio')+ xlab('')+
      annotate(geom="text", x='Down regulated', y=0.55, label='60',color="darkblue",size = 15)+
      annotate(geom="text", x='Down regulated', y=0.10, label='13',color="darkblue",size = 15)+
      annotate(geom="text", x='Up regulated', y=0.85, label='16',color="darkblue",size = 15)+
      annotate(geom="text", x='Up regulated', y=0.35, label='39',color="darkblue",size = 15)+
      lghplot.addthemeA(legend.position = 'top',size = 22,sizex = 20, sizey = 20)
dev.off()
```


60



<table>
<caption>A data.frame: 4 × 3</caption>
<thead>
	<tr><th scope=col>regulateType</th><th scope=col>metaboType</th><th scope=col>number</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Up regulated  </td><td>FAL    </td><td>16</td></tr>
	<tr><td>Up regulated  </td><td>non_FAL</td><td>39</td></tr>
	<tr><td>Down regulated</td><td>FAL    </td><td>60</td></tr>
	<tr><td>Down regulated</td><td>non_FAL</td><td>13</td></tr>
</tbody>
</table>




<strong>png:</strong> 2


### 5.5 DAscore analysis Figure 3C


```R
pdf('./figures/for_pub/Figure 3C.pdf',width = 12,height = 10)
lmMet$log2FC = lmMet$beta.cen
lmMet$Pvalue = lmMet$p.cen
lmMet$subSystemes = lmMet$metClass
DAscoreMet = Flux.subsystem_WDA_bootstrap(lmMet,bootstraplmMet[[1]],bootstraplmMet[[2]],lmMet,cutoff.nrxns = 3)
plot_DAscore(DAscoreMet$DAscore,subsystem = DAscoreMet$subsystem,fdr = DAscoreMet$fdr,
            gsize = 28,gsizex = 22,gsizey = 24)
dev.off()
```


<strong>png:</strong> 2


###  5.6 lipide changes in Cen Figure 3D


```R
# lipid data
lipidtype = c('PA','PC','PE','SM','PI')
mtype = substr(lmMet$met_name,1,2)
ida =  is.element(mtype,lipidtype) & lmMet$p.cen < 0.05
DElipid = lmMet[ida,]
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
pdf("./figures/for_pub/Figure 3D.pdf",width = 4,height = 8)

ggplot(tgc1, aes(x=met_name, y=log2fc)) + 
    geom_bar(aes(fill = tcolor),position=position_dodge(0.8),color = 'gray10',stat="identity",alpha = 0.7) + 
    geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.8,size = 0.8,alpha = 0.7,  
                  position=position_dodge(.9))+theme_bw() +coord_flip()+                  
lghplot.addtheme(hjust = F,size =12)+ xlab("") +  labs(y=expression(Log[2]("Fold change"))) +scale_fill_manual(values=c('darkblue','Brown'))
dev.off()

```


<ol class=list-inline>
	<li>71</li>
	<li>159</li>
</ol>




<strong>png:</strong> 2


### 5.7  Free fatty acid boxplot: Figure 3E


```R
pdf('./figures/for_pub/Figure 3E_trans_Vaccenic acid.pdf',width = 5,height = 6)
id = lmMet$met_name == 'trans-Vaccenic acid' & lmMet$score > 0.8
ida = clin.met$Class != 'F1'
p = ggplot(,aes(x = clin.met$Class[ida],y = met.expr[id,ida],color = clin.met$Class[ida],fill = clin.met$Class[ida]))+
lghplot.addthemeA(size = 24) + labs(y=expression(Log[2]("Peak Area")))  + xlab('')+labs(title = "Trans-Vaccenic")
print(lghplot.boxplot(p))
dev.off()
pdf('./figures/for_pub/Figure 3E_Palmitic acid.pdf',width = 5,height = 6)
id = lmMet$met_name == 'Palmitic acid' & lmMet$score > 0.8  & clin.met$Class != 'F1'
p = ggplot(,aes(x = clin.met$Class[ida],y = met.expr[id,ida],color = clin.met$Class[ida],fill = clin.met$Class[ida]))+
lghplot.addthemeA(size = 24) + ylab('') + xlab('')+labs(title = "Palmitic acid")
print(lghplot.boxplot(p))
dev.off()
```


<strong>png:</strong> 2



<strong>png:</strong> 2


### 5.8 clinical lipid Figure S6


```R
#### clin.met 
bx = lm(clin.met$TC ~ (clin.met$Class=='C') +  (clin.met$Gender == 'Female'),subset = clin.met$Class != 'F1')
tpval = car::Anova(bx)$`Pr(>F)`[1]
tbeta = bx$coefficients[2]

pdf("./figures/for_pub/Figure S6.pdf")
text = paste('p = ',signif(tpval,2),'\n','beta = ',signif(tbeta,2),sep = '')
id = clin.met$Class != 'F1'
p= ggplot(,aes(x = clin.met$Class[id],y = clin.met$TC[id], fill = clin.met$Class[id], color = clin.met$Class[id])) +
lghplot.addtheme()+ xlab('') + ylab('TC (umol/L)') + ggtitle('Serum total cholesterol') + 
annotate(geom="text", x='C', y=9, label=text,color="darkblue",size = 12)
print(lghplot.boxplot(p))
dev.off()
```


<strong>png:</strong> 2


##  5 Centenarians' plasma proteomics analysis


```R
protein = file2frame('./metabolic_modeling_centenarians/proteomics_270_matarix_simplifyA.txt')
protein = protein[,-c(2,3)]
protein.matrix = as.matrix(protein[,-1])
protein.genes = protein$Gene.name
colnames(protein.matrix) = sub('S','X',colnames(protein.matrix))
clin.protein = file2frame('./metabolic_modeling_centenarians/Hainan_CMWN_Clinical_simple.txt',row.names = 1)
list[IA,IB] = ismember(colnames(protein.matrix),rownames(clin.protein))
protein.matrix =protein.matrix[,IA]
clin.protein = clin.protein[IB,]
rownames(protein.matrix) = paste('protein_',1:nrow(protein.matrix),sep = '')
```


```R
DEprotein = DEGenes.simplified(protein.matrix,catagory = clin.protein$Class=='C',subset = clin.protein$Class != 'F1' & clin.protein$Gender== 'Female')
DEprotein$protein = protein.genes
```

### 5.1 boxplot ADIPOQ  Figure 4A


```R
id = protein.genes == 'ADIPOQ'
ADIPOQ = protein.matrix[id,]
#text = paste('p = ',signif(lmflux3['HMR_3406',]$p.cen,2),'\n','beta = ',signif(lmflux3['HMR_3406',]$beta.cen,2),sep = '')
list[ia,ib] = ismember(rownames(clin.protein),colnames(flux3.log2))
id = (clin.protein$Class == 'C' | clin.protein$Class == 'F1SP') & ia


bx = t.test(ADIPOQ~clin.protein$Class == 'C',subset = clin.protein$Class != 'F1' & id)
text = paste('p = ',signif(bx$p.value,2),sep = '')
pdf("./figures/for_pub/Figure 4A.pdf")
ggplot(,aes(clin.protein$Class[id],ADIPOQ[id],color = clin.protein$Class[id])) +geom_boxplot(outlier.size = -1) + 
      geom_jitter(height = 0,size =2) + 
      lghplot.addthemeA(sizex = 24,sizey = 22,size = 24) + ylim(0,0.4)+ xlab('')+ylab('ADIPOQ relative concentration')+
      annotate(geom="text", x='F1SP', y=0.3, label=text,color="darkblue",size = 10)
dev.off()
```


<strong>png:</strong> 2


### 5.2 lm model for age and gender effect: Figure 4B


```R
pdf("./figures/for_pub/Figure 4B.pdf")
id = protein.genes == 'ADIPOQ'
ADIPOQ = protein.matrix[id,]
list[ia,ib] = ismember(rownames(clin.protein),colnames(flux3.log2))
bx1 = summary(aov(ADIPOQ~clin.protein$Age + clin.protein$Gender, subset = clin.protein$Class != 'C' & ia))
text = paste('      Age: p = ',signif(bx1[[1]]$`Pr(>F)`[1],2),'\n','Gender: p = ',signif(bx1[[1]]$`Pr(>F)`[2],2),sep = '')
ida = clin.protein$Class!= 'C' & ia
ggplot(,aes(clin.protein$Age[ida],ADIPOQ[ida])) + geom_point(size = 2)+geom_smooth()+ ylim(0,0.4)+
      lghplot.addthemeA(sizex = 24,sizey = 22,size = 24)+ xlab('Age(years)')+ylab('ADIPOQ relative concentration')+
      annotate(geom="text", x=60, y=0.3, label=text,color="darkblue",size = 10)
dev.off()
```

    `geom_smooth()` using method = 'loess' and formula 'y ~ x'
    


<strong>png:</strong> 2


### 6.3 Figure S7

#### Figure S7A


```R
agmatine = file2frame('./metabolic_modeling_centenarians/agmatine_from_GC.txt')
rownames(agmatine) = paste0('X',agmatine$ID)
clin.agmatine = file2frame('./metabolic_modeling_centenarians/Hainan_CMWN_Clinical_simple.txt',row.names = 1)
clin.agmatine = clin.agmatine[rownames(agmatine),]
pdf("./figures/for_pub/Figure S7A.pdf")
id = (clin.agmatine$Class == 'C' | clin.agmatine$Class == 'F1SP')

bx = t.test(agmatine$agmatine~clin.agmatine$Class == 'C',subset = clin.agmatine$Class != 'F1' & id)
text = paste('p = ',signif(bx$p.value,2),sep = '')

p = ggplot(,aes(clin.agmatine$Class[id],agmatine$agmatine[id],color = clin.agmatine$Class[id])) + 
      lghplot.addthemeA(sizex = 24,sizey = 22,size = 24) + xlab('')+ylab('Agmatine relative concentration')+
      annotate(geom="text", x='C', y=0.007, label=text,color="darkblue",size = 10)
print(lghplot.boxplot(p))
dev.off()
```


<strong>png:</strong> 2


#### Figure S7B


```R
pdf("./figures/for_pub/Figure S7B_1.pdf")
id = protein.genes == 'LEP'
LEP = protein.matrix[id,]
#text = paste('p = ',signif(lmflux3['HMR_3406',]$p.cen,2),'\n','beta = ',signif(lmflux3['HMR_3406',]$beta.cen,2),sep = '')
list[ia,ib] = ismember(rownames(clin.protein),colnames(flux3.log2))
id = (clin.protein$Class == 'C' | clin.protein$Class == 'F1SP') & ia
bx = t.test(LEP~clin.protein$Class == 'C',subset = clin.protein$Class != 'F1' & id)
text = paste('p = ',signif(bx$p.value,2),sep = '')

p = ggplot(,aes(clin.protein$Class[id],LEP[id],color = clin.protein$Class[id])) + 
      lghplot.addthemeA(sizex = 24,sizey = 22,size = 24) + ylim(0,0.4)+ xlab('')+ylab('LEP relative concentration')+
      annotate(geom="text", x='C', y=0.35, label=text,color="darkblue",size = 10)
print(lghplot.boxplot(p))
dev.off()
```


<strong>png:</strong> 2

