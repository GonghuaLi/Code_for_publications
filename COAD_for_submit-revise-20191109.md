
## 1 set pwd and subroutines


```R
options(warn = -1)    
setwd("D:/Work/coll/申秋/GDE1-revise")
library(ggplot2)
library(survival)
library(survminer)
library(Rmisc)
library(reshape2)
```


```R
#subroutines
file2frame <- function(file, header = TRUE, sep = "\t",stringsAsFactors = F,row.names = NULL){
    out = read.delim(file, header = header, sep = sep,stringsAsFactors = stringsAsFactors,row.names = row.names)
    return(out)
}
as.numeric.matrix <- function(s){
  v = apply(s,2,as.numeric)
  rownames(v) = rownames(s)
  colnames(v) = colnames(s)
  return(v)
}

list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
    args <- as.list(match.call())
    args <- args[-c(1:2,length(args))]
    length(value) <- length(args)
    for(i in seq(along=args)) {
        a <- args[[i]]
        if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
    }
    x
}

ismember <- function(A,B){
    IA = is.element(A,B)
    IB = matrix(0,length(A),1)
    for (i in 1:length(A)) {
        if (IA[i]){
            c= which(B==A[i]);
            IB[i] = c[1]
        }
    }
    IB = as.vector(IB)
    return(list(IA,IB))
}

#plot
lghplot.boxplot <- function(p,width = 0.3){
  p = p +geom_violin(alpha = 0.5,outlier.size = 0) +
    geom_boxplot(width = width,color = 'black',fill = 'white',outlier.size = 0,alpha = 1)

}
lghplot.addtheme <- function (size = 18,hjust = FALSE,legend.position = "none"){
 p = theme(axis.text.x = element_text(size = size,face = 'bold'))+
    theme(axis.text.y = element_text(size = size,face = 'bold'))+
    theme(axis.title=element_text(size=size,face = 'bold')) +
    theme(axis.title=element_text(size=size,face = 'bold')) +
    theme(title = element_text(size= size,face = 'bold'))+
    theme(#panel.grid.major =element_blank(),
        strip.text = element_text(size=size,face = 'bold'),
        legend.position  = legend.position,
       #legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
          )
 if (hjust){
   p = p +theme(axis.text.x = element_text(angle = 45, hjust = 1))
 }
 return(p)
}

log2fc <- function(e,c){
    return(mean(e[c])-mean(e[!c]))
}
pval <- function(e,c){
    xx = summary(aov(e~c))
    return(xx[[1]]$`Pr(>F)`[1])
}
```

## 2 read data


```R
#readnorm
expr.norm = file2frame('./data/expr_norm_tpm.txt')
expr.norm = expr.norm[!duplicated(expr.norm$gene),]
rownames(expr.norm) = expr.norm$gene
expr.norm = as.numeric.matrix(expr.norm[,-1])

#read case
expr.case = file2frame('./data/expr_case_tpm.txt')
expr.case = expr.case[!duplicated(expr.case$gene),]
rownames(expr.case) = expr.case$gene
expr.case = as.numeric.matrix(expr.case[,-1])

#clin
clin = file2frame('./data/coad_clin_data.csv',row.names =1)
clin = clin[clin$HISTOLOGICAL_DIAGNOSIS != '',]
list[IA,IB] = ismember(rownames(clin),colnames(expr.case))
clin = clin[IA,]
expr.case = expr.case[,IB]

expr.all = cbind(expr.norm,expr.case)
```

## 2 Sample QC


```R
sampleclass = c(rep('green',ncol(expr.norm)),rep('pink',ncol(expr.case)))
pdf('./samples_boxplot1.pdf')
boxplot(expr.all,outline = F,outwex = 0,xlab = 'Sample ID',ylab = 'log2 TPM',names = 1:ncol(expr.all),col = sampleclass,
        pars = list(boxwex = 1, staplewex =1, outwex =0.5),
          horizontal = FALSE, add = FALSE, at = NULL)
dev.off()
```


<strong>png:</strong> 2


## 3 Figures

### 3.1 Figure 4A


```R
sum(tclass_simple == 'Normal')
sum(tclass_simple == 'N-COAD')
sum(tclass_simple == 'M-COAD')
```


41



233



38



```R
tclass = c(rep('Control',ncol(expr.norm)),clin$HISTOLOGICAL_DIAGNOSIS)
tclass_simple = tclass
tclass_simple[tclass_simple == 'Colon Adenocarcinoma'] = 'N-COAD'
tclass_simple[tclass_simple == 'Colon Mucinous Adenocarcinoma'] = 'M-COAD'
tclass_simple[tclass_simple == 'Control'] = 'Normal'
tclass_simple = factor(tclass_simple,levels = c('N-COAD','M-COAD','Normal'))
idaa = 'GDE1'
#pdf('./figures/Figure4A_GDE1_boxplot.pdf')
p = ggplot(,aes(x = tclass_simple,y = expr.all[idaa,],color = tclass_simple))+ ylab('log2(TPM)')+
    geom_text(x=3, y=30, label="Scatter plot")
    xlab('')+theme_bw()+lghplot.addtheme(size = 22)
print(lghplot.boxplot(p))
#dev.off()

```


    NULL



![png](output_10_1.png)


### 3.2 Figure 4B


```R
#COAD vs COMAD
tgene = 'GDE1'
ida = clin$HISTOLOGICAL_DIAGNOSIS != ''
cutoff = median(expr.case[tgene,ida])
survival_time = clin$OS_MONTHS
survival_status = clin$OS_STATUS == 'DECEASED'
st0 = clin$HISTOLOGICAL_DIAGNOSIS == 'Colon Adenocarcinoma'
st = rep('M-COAD',length(st0))
st[st0] = 'N-COAD'
tdata = data.frame(cancer_type = st[ida],
                  survival_time = survival_time[ida],
                  survival_status = survival_status[ida])
survdiff(Surv(survival_time[ida],survival_status[ida])~st[ida])

tfit = survfit(Surv(survival_time,survival_status)~cancer_type,data = tdata)
p = ggsurvplot(tfit,legend.labs = c("M-COAD", "N-COAD"),
           risk.table.col = "strata", 
           ggtheme = theme_bw(base_size = 22), 
           palette = c("#CC0000", "#4C9900"),
           xlim = c(0, 120))
#pdf('./figures/Figure4B_GDE1_N-COADvsM-COAD_survival.pdf')
print(p)
#dev.off()
```


    Call:
    survdiff(formula = Surv(survival_time[ida], survival_status[ida]) ~ 
        st[ida])
    
                     N Observed Expected (O-E)^2/E (O-E)^2/V
    st[ida]=M-COAD  38        9     8.02    0.1197     0.138
    st[ida]=N-COAD 233       53    53.98    0.0178     0.138
    
     Chisq= 0.1  on 1 degrees of freedom, p= 0.7 



![png](output_12_1.png)


### 3.3 Figure 4C


```R
# N-COAD
tgene = 'GDE1'
ida = clin$HISTOLOGICAL_DIAGNOSIS == 'Colon Adenocarcinoma'
cutoff = median(expr.case[tgene,ida])
survival_time = clin$OS_MONTHS
survival_status = clin$OS_STATUS == 'DECEASED'
st0 = expr.case[tgene,] > cutoff
st = rep('low expression',length(st0))
st[st0] = 'high expression'
tdata = data.frame(cancer_type = st[ida],
                  survival_time = survival_time[ida],
                  survival_status = survival_status[ida])
survdiff(Surv(survival_time[ida],survival_status[ida])~st[ida])

tfit = survfit(Surv(survival_time,survival_status)~cancer_type,data = tdata)
p = ggsurvplot(tfit,legend.labs = c("high expression", "low expression"),
           risk.table.col = "strata", 
           ggtheme = theme_bw(base_size = 22), 
           palette = c("#CC0000", "#4C9900"),
           xlim = c(0, 120))
#pdf('./figures/Figure4B_GDE1_N-COAD_survival.pdf')
print(p)
#dev.off()
```


    Call:
    survdiff(formula = Surv(survival_time[ida], survival_status[ida]) ~ 
        st[ida])
    
                              N Observed Expected (O-E)^2/E (O-E)^2/V
    st[ida]=high expression 116       15     26.4      4.90      9.86
    st[ida]=low expression  117       38     26.6      4.85      9.86
    
     Chisq= 9.9  on 1 degrees of freedom, p= 0.002 



![png](output_14_1.png)


### 3.4 Figure 4D


```R
#COMAD
tgene = 'GDE1'
ida = clin$HISTOLOGICAL_DIAGNOSIS != 'Colon Adenocarcinoma'
cutoff = median(expr.case[tgene,ida])
survival_time = clin$OS_MONTHS
survival_status = clin$OS_STATUS == 'DECEASED'
st0 = expr.case[tgene,] > cutoff
st = rep('low expression',length(st0))
st[st0] = 'high expression'
tdata = data.frame(cancer_type = st[ida],
                  survival_time = survival_time[ida],
                  survival_status = survival_status[ida])
survdiff(Surv(survival_time[ida],survival_status[ida])~st[ida])
tfit = survfit(Surv(survival_time,survival_status)~cancer_type,data = tdata)
p = ggsurvplot(tfit,legend.labs = c("high expression", "low expression"),
           risk.table.col = "strata",
           ggtheme = theme_bw(base_size = 22),
           palette = c("#CC0000", "#4C9900"),
           xlim = c(0, 120))
#pdf('./figures/Figure4B_GDE1_M-COMD_survival.pdf')
print(p)
#dev.off()
```


    Call:
    survdiff(formula = Surv(survival_time[ida], survival_status[ida]) ~ 
        st[ida])
    
                             N Observed Expected (O-E)^2/E (O-E)^2/V
    st[ida]=high expression 19        5     5.15   0.00458    0.0111
    st[ida]=low expression  19        4     3.85   0.00614    0.0111
    
     Chisq= 0  on 1 degrees of freedom, p= 0.9 



![png](output_16_1.png)


## 4 Tables


```R
des = file2frame('./data/des.all_genes.txt',row.names = 1) # supplementary figure 1
ogenes = file2frame('./data/WGCNAmodule_revise_red.txt',row.names = NULL)
rownames(ogenes) = ogenes$gene_name
ogenes$log2fc = rep(1,nrow(ogenes))
ogenes$pvalue = rep(1,nrow(ogenes))
ogenes$p.adjust = rep(1,nrow(ogenes))
ogenes$HR_All_COAD =rep(1,nrow(ogenes))
ogenes$survival.pAll_COAD = rep(1,nrow(ogenes))
ogenes$out_All_COAD = rep('x',nrow(ogenes))
ogenes$HR_N_COAD =rep(1,nrow(ogenes))
ogenes$survival.pN_COAD = rep(1,nrow(ogenes))
ogenes$out_N_COAD = rep('x',nrow(ogenes))
ogenes$HR_M_COAD =rep(1,nrow(ogenes))
ogenes$survival.pM_COAD = rep(1,nrow(ogenes))
ogenes$out_M_COAD = rep('x',nrow(ogenes))

for (i in 1:nrow(ogenes)){
    
    # log2fc
    tgene = ogenes$gene_name[i]
    
    # all
    ida = clin$HISTOLOGICAL_DIAGNOSIS != ' '
    cutoff = median(expr.case[tgene,ida])
    survival_time = clin$OS_MONTHS
    survival_status = clin$OS_STATUS == 'DECEASED'
    st = expr.case[tgene,] > cutoff
    xx = survdiff(Surv(survival_time[ida],survival_status[ida])~st[ida])
    ogenes$survival.pAll_COAD[i] = pchisq(xx$chisq, length(xx$n)-1, lower.tail = FALSE)
    ogenes$HR_All_COAD[i] = (xx$obs[2]/xx$exp[2])/(xx$obs[1]/xx$exp[1])
    ogenes$out_All_COAD[i] = paste(as.character(signif(ogenes$HR_All_COAD[i],2)),'(',
                            as.character(signif(ogenes$survival.pAll_COAD[i],2)),')',sep='')
    
    #N-COAD
    ida = clin$HISTOLOGICAL_DIAGNOSIS == 'Colon Adenocarcinoma'
    cutoff = median(expr.case[tgene,ida])
    survival_time = clin$OS_MONTHS
    survival_status = clin$OS_STATUS == 'DECEASED'
    st = expr.case[tgene,] > cutoff
    xx = survdiff(Surv(survival_time[ida],survival_status[ida])~st[ida])
    ogenes$survival.pN_COAD[i] = pchisq(xx$chisq, length(xx$n)-1, lower.tail = FALSE)
    ogenes$HR_N_COAD[i] = (xx$obs[2]/xx$exp[2])/(xx$obs[1]/xx$exp[1])
    ogenes$out_N_COAD[i] = paste(as.character(signif(ogenes$HR_N_COAD[i],2)),'(',
                            as.character(signif(ogenes$survival.pN_COAD[i],2)),')',sep='')
    
    #M-COAD
    ida = clin$HISTOLOGICAL_DIAGNOSIS != 'Colon Adenocarcinoma'
    cutoff = median(expr.case[tgene,ida])
    survival_time = clin$OS_MONTHS
    survival_status = clin$OS_STATUS == 'DECEASED'
    st = expr.case[tgene,] > cutoff
    xx = survdiff(Surv(survival_time[ida],survival_status[ida])~st[ida])
    ogenes$survival.pM_COAD[i] = pchisq(xx$chisq, length(xx$n)-1, lower.tail = FALSE)
    ogenes$HR_M_COAD[i] = (xx$obs[2]/xx$exp[2])/(xx$obs[1]/xx$exp[1])
    ogenes$out_M_COAD[i] = paste(as.character(signif(ogenes$HR_M_COAD[i],2)),'(',
                            as.character(signif(ogenes$survival.pM_COAD[i],2)),')',sep='')
    
}
ix = sort.int(ogenes$survival.pAll_COAD,decreasing = F,index.return = T)$ix
ogenes = ogenes[ix,]
ogenes$log2fc = des[rownames(ogenes),]$log2FC
ogenes$pvalue = des[rownames(ogenes),]$pvalue
ogenes$p.adjust = des[rownames(ogenes),]$p.adjust
write.table(ogenes,'Table 1_revise.txt',sep = '\t',row.names = F,col.names = T,quote = F)
```


```R
head(ogenes)
```


<table>
<caption>A data.frame: 6 × 14</caption>
<thead>
	<tr><th></th><th scope=col>gene_name</th><th scope=col>modeduleName</th><th scope=col>log2fc</th><th scope=col>pvalue</th><th scope=col>p.adjust</th><th scope=col>HR_All_COAD</th><th scope=col>survival.pAll_COAD</th><th scope=col>out_All_COAD</th><th scope=col>HR_N_COAD</th><th scope=col>survival.pN_COAD</th><th scope=col>out_N_COAD</th><th scope=col>HR_M_COAD</th><th scope=col>survival.pM_COAD</th><th scope=col>out_M_COAD</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>GDE1</th><td>GDE1   </td><td>red</td><td>-0.5787162</td><td>3.29e-21</td><td>2.10e-20</td><td>0.4118909</td><td>0.0009868335</td><td>0.41(0.00099)</td><td>0.3986529</td><td>0.001687905</td><td>0.4(0.0017)</td><td>0.9329231</td><td>0.91608648</td><td>0.93(0.92)</td></tr>
	<tr><th scope=row>DARS</th><td>DARS   </td><td>red</td><td> 1.1386487</td><td>4.40e-52</td><td>2.63e-50</td><td>2.0597543</td><td>0.0043170964</td><td>2.1(0.0043)  </td><td>2.0076648</td><td>0.009103973</td><td>2(0.0091)  </td><td>2.9458436</td><td>0.09954682</td><td>2.9(0.1)  </td></tr>
	<tr><th scope=row>SLCO4A1</th><td>SLCO4A1</td><td>red</td><td> 3.4952825</td><td>5.42e-60</td><td>8.27e-58</td><td>1.9828334</td><td>0.0065671052</td><td>2(0.0066)    </td><td>1.7417607</td><td>0.042882423</td><td>1.7(0.043) </td><td>3.0917807</td><td>0.08674171</td><td>3.1(0.087)</td></tr>
	<tr><th scope=row>EPHB2</th><td>EPHB2  </td><td>red</td><td> 1.1015325</td><td>2.44e-24</td><td>1.84e-23</td><td>0.5341220</td><td>0.0162134563</td><td>0.53(0.016)  </td><td>0.5103084</td><td>0.019434768</td><td>0.51(0.019)</td><td>0.9669079</td><td>0.95743045</td><td>0.97(0.96)</td></tr>
	<tr><th scope=row>B3GNT8</th><td>B3GNT8 </td><td>red</td><td>-1.6851371</td><td>3.15e-13</td><td>1.02e-12</td><td>0.5449521</td><td>0.0199357953</td><td>0.54(0.02)   </td><td>0.5250038</td><td>0.021917458</td><td>0.53(0.022)</td><td>0.6814847</td><td>0.56013366</td><td>0.68(0.56)</td></tr>
	<tr><th scope=row>DGUOK</th><td>DGUOK  </td><td>red</td><td> 0.3719217</td><td>1.92e-14</td><td>6.96e-14</td><td>1.7411089</td><td>0.0292847406</td><td>1.7(0.029)   </td><td>1.7123068</td><td>0.051907061</td><td>1.7(0.052) </td><td>2.2587654</td><td>0.19976475</td><td>2.3(0.2)  </td></tr>
</tbody>
</table>


