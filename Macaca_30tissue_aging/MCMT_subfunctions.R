standardise_matrix <- function (m) 
{
    out = m
    for (i in 1:nrow(m)) {
        out[i, ] <- (m[i, ] - mean(m[i, ], na.rm = TRUE))/sd(m[i,], na.rm = TRUE)
    }
    out = out[-1,]
    return(out)
}

standardise_matrix_1 <- function (m) 
{
    out = m
    for (i in 1:nrow(m)) {
        out[i, ] <- (m[i, ] - m[i, 1])/sd(m[i,], na.rm = TRUE)
    }
    out = out[-1,]
    return(out)
}

rbind2 <- function(m1, m2){
    m1f = as.data.frame(m1)
    #rownames(m1f) = rownames(m1f)
    m2f = as.data.frame(m2)
    #rownames(m2f) = rownames(m2f)
    out = as.matrix(as.data.frame(rbindlist(list(m1f,m2f),fill = T)))
    rownames(out)  =c(rownames(m1f),rownames(m2f))
    return(out)
}


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}

#
shorten_names_forGO1 <- function (x, n_word = 5, n_char = 30) 
{
    if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 
        n_char)) {
        if (nchar(x) > n_char) 
            x <- substr(x, 1, n_char)
        x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x, 
            " ")[[1]]), n_word)], collapse = " "), "...", sep = "")
        return(x)
    }
    else {
        return(x)
    }
}
plot_enrich <- function(goup, godown,n_char = 60){
    CPCOLS =c('#BB0021','#3B4992')
    aa = sort(godown[1:min(10,length(godown))],decreasing = F)
    bb = sort(goup[1:min(10,length(goup))],decreasing = F)
    tgoname = capitalize(c(names(aa),names(bb)))
    tgoname[duplicated(tgoname)] = paste0(' ',tgoname[duplicated(tgoname)])
    tgoname.short = as.vector(sapply(tgoname, shorten_names_forGO1,n_char = n_char))
    tmp = data.frame(log10pval = as.vector(c(aa,bb)),
                   direction = factor(c(rep('Down',length(aa)),rep('Up',length(bb))),
                                      levels = c('Up','Down')),
                  goname =   factor(tgoname, levels = tgoname),
                    showname = factor(tgoname.short,levels = tgoname.short)
                    )
   # print(tmp)
    #tlab = (sapply(levels(tmp$goname), shorten_names_forGO))
    tlab = tmp$showname
    tmp$number <- factor(1:nrow(tmp))
     p = ggplot(data = tmp, aes(x = number, y = log10pval, fill = direction)) + 
         geom_bar(stat = "identity", width = 0.8) + coord_flip() +
         scale_fill_manual(values = CPCOLS) + theme_classic() + 
        scale_x_discrete(labels = tlab) + xlab("GO term") + theme(axis.text = element_text(face = "bold", 
        color = "gray0")) + theme(axis.text.y = element_text(size = 8,face = "bold"))+lghplot.addtheme(size = 10)+
        ylab('-log10(p.adjust)')+ xlab('')
    return(p)
    
}

fmt_dcimals <- function(decimals=0){
    function(x) format(x,nsmall = decimals,scientific = FALSE)
}

# sub routine plot volcano
plot_Volcano <- function(res2, title){
    

  tmpup = res2$Pvalue
  tmpup[is.na(tmpup)] = 1
  tmpup[res2$log2FC < 0] = 1
  sortid_up = sort.int(-log10(tmpup),decreasing = T,index.return = T)$ix
  tmpid.up = res2$ID[sortid_up[1:5]]

  tmpdown = res2$Pvalue
  tmpdown[is.na(tmpdown)] = 1
  tmpdown[res2$log2FC > 0] = 1
  sortid_down = sort.int(-log10(tmpdown),decreasing = T,index.return = T)$ix
  tmpid.down = res2$ID[sortid_down[1:5]]

  vid = is.element(res2$ID,c(tmpid.up,tmpid.down))
  tlab = res2$ID
  tlab[!vid] = NA

  keyvals <- rep('gray50', nrow(res2))
  names(keyvals) <- rep('NS', nrow(res2))   

  keyvals[which(res2$log2FC > 0.58 & res2$Pvalue < 0.05)] <- "Brown"    
  names(keyvals)[which(res2$log2FC > 0.58 & res2$Pvalue < 0.05)] <- 'High'    

  keyvals[which(res2$log2FC < -0.58 & res2$Pvalue < 0.05)] <- "darkblue"    
  names(keyvals)[which(res2$log2FC < -0.58 & res2$Pvalue < 0.05)] <- 'Low'   
  p = EnhancedVolcano(res2,
                        lab = tlab,
                        x = 'log2FC',
                        y = 'Pvalue',
                        caption = NULL,
                        title = title,
                        border = 'full',
                        titleLabSize = 12,
                        FCcutoff = 0.58,
                        cutoffLineWidth = F,
                        axisLabSize = 12,
                        subtitle = NULL,
                        cutoffLineCol = 'white',
                        gridlines.minor = F,
                        gridlines.major = F,
                        xlab = bquote(~Log[2]~ 'fold change'),
                        pCutoff = 0.05,
                        colCustom = keyvals,
                        colAlpha = 4/5,
                        legendPosition = 'none',
                        legendLabSize = 5,
                        legendIconSize = 3,
                        drawConnectors = TRUE,
                        widthConnectors = 0.5,
                        pointSize = -0.3*log10(res2$Pvalue),labSize = 3,
                        colConnectors = 'black')
    return(p)
}


list_to_matrix <- function(DEproFC,alltissues){
    DEproFC_matrix = list()
    for(i in 1:length(alltissues)){
        tmp = matrix(DEproFC[[i]],1,length(DEproFC[[i]]))
        tmp = as.data.frame(tmp)
        colnames(tmp) = names(DEproFC[[i]])
        DEproFC_matrix[[i]] = tmp
    }
    DEproFC_matrix = t(as.matrix(rbindlist(DEproFC_matrix,fill = T)))
    colnames(DEproFC_matrix) = names(DEproFC)
    #vid = rowSums(is.na(DEproFC_matrix)) < ncol(DEproFC_matrix)/2
    #DEproFC_matrix = DEproFC_matrix[vid,]
    return(DEproFC_matrix)
}

met.class_enrichment <- function(mets,annote){
  require(clusterProfiler)

  vmet = intersect(mets,rownames(annote))
  fluxgmt = data.frame(ont = annote$sub_class,
                       gene = rownames(annote),stringsAsFactors = F) 
  
  Recon3D <- enricher(gene = vmet,
                TERM2GENE=fluxgmt,
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =1,
                qvalueCutoff = 1
                )
  Recon3Dout = Recon3D@result
  Recon3Dout$DB = rep('HMDBclass',dim(Recon3Dout)[1])
  n = dim(Recon3Dout)[2]
  Recon3Dout = Recon3Dout[,c(n,1:(n-1))]
  return(Recon3Dout)
}

fillgaps_rowMinNA_1 <- function(X){
  out = X
  for(i in 1:nrow(out)){
    id = which(is.na(out[i,]))
    if(!any(id)){next;}
    e = min(out[i,],na.rm = T)
    out[i,id] = e
  }
  out[is.infinite(out)] = NA
  return(out)
}

fillgaps_MinNA <- function(X){
  out = X
  min_e = min(out,na.rm =T)
  alpha = max(min_e,1)
  for(i in 1:nrow(out)){
    id = which(is.na(out[i,]))
    if(!any(id)){next;}
    #e = min(out[i,],na.rm = T)
    #if(is.infinite(e)){
    out[i,id] = min_e  + alpha*runif(n=length(id), min = -0.01, max=0.01)
    #}else{
    #  out[i,id] = e
    #}
    
  }
  return(out)
}

fillgaps_MinNA_formeta <- function(X,nas,nascutoff = 0.3){
  out = X
  min_e = min(out,na.rm =T)
  #alpha = max(min_e,1)
  for(i in 1:nrow(out)){
    id = which(is.na(out[i,]))
    if(!any(id)){next;}
    #e = min(out[i,],na.rm = T)
    if(nas[i] < nascutoff){
      #e = min(out[i,],na.rm = T)
      out[i,id] = min_e*(1+runif(n=length(id), min = -0.01, max=0.01))
    }
  }
  return(out)
}


lmGene <- function(repExpr,repClin){
  # set data
  flux3.log2 = list()
  
  allrxns = c()
  thenames = names(repClin)
  
  
  for(i in 1:length(thenames)){
    studyname = thenames[i]
    #expr.log2 = #log2(repExpr[[studyname]]+1e-6)
    expr.log2 = repExpr[[studyname]]
    allrxns = unique(c(allrxns,rownames(expr.log2)))
    flux3.log2[[i]] = expr.log2
  }
  # filter allrxns NA
  nas  = rep(0,length(allrxns))
  k = 0
  for(i in 1:length(thenames)){
    tmp = as.data.frame(flux3.log2[[i]])
    tmp = as.matrix(tmp[allrxns,])
    k = k+ncol(tmp)
    nas = nas + rowSums(is.na(tmp))
  }
  nas = nas/k
  names(nas) = allrxns
  
  #allrxns = allrxns[nas/k < 0.5]
  
  flux3.AveExpr = matrix(0,length(allrxns),length(repClin))
  flux3.log2FC = matrix(0,length(allrxns),length(repClin))
  flux3.Pvalue = matrix(0,length(allrxns),length(repClin))
  flux3.FDR = matrix(0,length(allrxns),length(repClin))
  xtissues = c()
  xage  =c()
  
  for(i in 1:length(thenames)){
    tmp = as.data.frame(flux3.log2[[i]])
    tmp = as.matrix(tmp[allrxns,])
    tmp1 = fillgaps_rowMinNA_1(tmp)
    tmp = fillgaps_MinNA_formeta(tmp,nas)
    #tmp[is.na(tmp)] = 0
    rownames(tmp) = allrxns
    rownames(tmp1) = allrxns
    
    xtissues = c(xtissues,rep(thenames[i],ncol(tmp)))
    xage = c(xage,repClin[[i]]$age)
    
    tmpDEflux = DEGenes.simplified(tmp1,catagory = repClin[[i]]$type =='Elderly',
                                   subset = repClin[[i]]$type =='Elderly' | repClin[[i]]$type =='Young_adult')
    #tmp_filled = fillgaps_rowMinNA(tmp)
    #ss = norm::prelim.norm(tmp)
    #tmpfilled = norm::imp.norm(ss,'em.norm',tmp)
    #tmp_filled = fillgaps_rowMinNA(tmp)
    if(i ==1){
      exprall = tmp
    }else{
      exprall = cbind(exprall,tmp)
    }
      
    
    flux3.log2FC[,i] = tmpDEflux$log2FC
    flux3.Pvalue[,i] = tmpDEflux$Pvalue
    flux3.FDR[,i] = tmpDEflux$FDR
    flux3.AveExpr[,i]  =  tmpDEflux$AveExpr
  }
  
  
  colnames(flux3.AveExpr) = paste0('AveExpr_',thenames)
  colnames(flux3.log2FC) = paste0('log2FC_',thenames)
  colnames(flux3.Pvalue) = paste0('Pvalue_',thenames)
  colnames(flux3.FDR) = paste0('FDR_',thenames)
  
  lmflux3 = data.frame(fdr.age = rep(1,length(allrxns)),
                      p.age = rep(1,length(allrxns)),
                      beta.age = rep(0,length(allrxns)))

  
  exprall = as.matrix(exprall)
  for(i in 1:length(allrxns)){
    tmpx = exprall[i,]
    
    if(any(is.na(tmpx))){next;}
    
    bx = lm(tmpx ~ xage+ xtissues)
    lmflux3$beta.age[i] = bx$coefficients[2]
    lmflux3$p.age[i] = car::Anova(bx)$`Pr(>F)`[1]
  }
  
  lmflux3$fdr.age[lmflux3$p.age < 1] = p.adjust(lmflux3$p.age[lmflux3$p.age < 1],method = "fdr")
  
  #allrxns = allrxns[nas/k < 0.5]
  out = data.frame(ID = allrxns,stringsAsFactors = F,
                   MetaAveExpr = rowMeans(flux3.AveExpr,na.rm = T),
                   Metalog2FC = rowMeans(flux3.log2FC,na.rm = T),
                   p.age = lmflux3$p.age,
                   fdr.age = lmflux3$fdr.age,
                   beta.age = lmflux3$beta.age)
  
  out = cbind(out,flux3.log2FC)
  out = cbind(out,flux3.Pvalue)
  out = cbind(out,flux3.FDR)
  #flux3.annote = Recon3.annote
  #out = cbind(out,flux3.annote[allrxns,])
  rownames(out) = allrxns
  
  return(out)

  
}



metaGene <- function(repExpr,repClin, 
                     meta.method = "FEM",tail = 'abs',parametric = TRUE){
  # set data
  flux3.log2 = list()

  allrxns = c()
  thenames = names(repClin)
  
  for(i in 1:length(thenames)){
    studyname = thenames[i]
    #expr.log2 = #log2(repExpr[[studyname]]+1e-6)
    expr.log2 = repExpr[[studyname]]
    allrxns = unique(c(allrxns,rownames(expr.log2)))
    flux3.log2[[i]] = expr.log2
  }
  # filter allrxns NA
  nas  = rep(0,length(allrxns))
  k = 0
  for(i in 1:length(thenames)){
    tmp = as.data.frame(flux3.log2[[i]])
    tmp = as.matrix(tmp[allrxns,])
    k = k+ncol(tmp)
    nas = nas + rowSums(is.na(tmp))
  }
  nas = nas/k
  names(nas) = allrxns
  
  #allrxns = allrxns[nas/k < 0.5]
  
  flux3.AveExpr = matrix(0,length(allrxns),length(repClin))
  flux3.log2FC = matrix(0,length(allrxns),length(repClin))
  flux3.Pvalue = matrix(0,length(allrxns),length(repClin))
  flux3.FDR = matrix(0,length(allrxns),length(repClin))
  for(i in 1:length(thenames)){
    tmp = as.data.frame(flux3.log2[[i]])
    tmp = as.matrix(tmp[allrxns,])
    tmp1 = fillgaps_rowMinNA_1(tmp)
    tmp = fillgaps_MinNA_formeta(tmp,nas)
    #tmp[is.na(tmp)] = 0
    rownames(tmp) = allrxns
    rownames(tmp1) = allrxns
   
    tmpDEflux = DEGenes.simplified(tmp1,catagory = repClin[[i]]$type =='Elderly',
                                   subset = repClin[[i]]$type =='Elderly' | repClin[[i]]$type =='Young_adult')
    #tmp_filled = fillgaps_rowMinNA(tmp)
    #ss = norm::prelim.norm(tmp)
    #tmpfilled = norm::imp.norm(ss,'em.norm',tmp)
    #tmp_filled = fillgaps_rowMinNA(tmp)
    flux3.log2[[i]] = tmp
    flux3.log2FC[,i] = tmpDEflux$log2FC
    flux3.Pvalue[,i] = tmpDEflux$Pvalue
    flux3.FDR[,i] = tmpDEflux$FDR
    flux3.AveExpr[,i]  =  tmpDEflux$AveExpr
  }
  colnames(flux3.AveExpr) = paste0('AveExpr_',thenames)
  colnames(flux3.log2FC) = paste0('log2FC_',thenames)
  colnames(flux3.Pvalue) = paste0('Pvalue_',thenames)
  colnames(flux3.FDR) = paste0('FDR_',thenames)
  
  #run metaDE
  data = flux3.log2
  label = list()
  K <- length(data)
  for(i in 1:length(repClin)){
    label[[i]] = repClin[[i]]$type
  }
  clin.data <- lapply(label, function(x) {data.frame(x)} )
  for (k in 1:length(clin.data)){
  colnames(clin.data[[k]]) <- "label"
  }
  select.group <- c('Elderly','Young_adult')
  ref.level <- "Young_adult"
  data.type <- "continuous"
  ind.method <- rep('limma',length(data))
  resp.type <- "twoclass"
  paired <- rep(FALSE,length(data))
  meta.res <- MetaDE(data=data,clin.data = clin.data,
                    data.type=data.type,resp.type = resp.type,
                    response='label',
                    ind.method=ind.method, meta.method=meta.method,
                    select.group = select.group, ref.level=ref.level,
                    REM.type = 'HS',
                    paired=paired,tail=tail,parametric=parametric)
  #
  out = data.frame(ID = allrxns,stringsAsFactors = F,
                   MetaAveExpr = rowMeans(flux3.AveExpr,na.rm = T),
                   Metalog2FC = rowMeans(flux3.log2FC,na.rm = T),
                   MetaPvalue = as.vector(meta.res$meta.analysis$pval),
                   MetaFDR = as.vector(meta.res$meta.analysis$FDR),
                   MetaZvalue = as.vector(meta.res$meta.analysis$zval))
  out$MetaFDR[is.na(out$MetaFDR)] = 1
  out$MetaPvalue[is.na(out$MetaPvalue)] = 1
  out$MetaZvalue[is.na(out$MetaZvalue)] = 0
  
  out = cbind(out,flux3.log2FC)
  out = cbind(out,flux3.Pvalue)
  out = cbind(out,flux3.FDR)
  #flux3.annote = Recon3.annote
  #out = cbind(out,flux3.annote[allrxns,])

  # re-calculate pvalue and FDR for pvalue < eps
  idx = out$MetaPvalue < 1.0
  tmpp = pnorm(-abs(out$MetaZvalue[idx]))*2
  tmpfdr = p.adjust(tmpp,method = "BH")
  out$MetaFDR[idx] = tmpfdr
  out$MetaPvalue[idx] = tmpp
  
  rownames(out) = allrxns
  return(out)
  
}

list_element_select <- function(x, rname,cname){
  out = matrix('',length(rname),length(x))
  for(i in 1:length(x)){
    out[,i] = unlist(x[[i]][rname,cname])
  }
  colnames(out) = names(x)
  rownames(out) = rname
  return(out)
}


get_corr <- function(list1,list2){
  outcor = list()
  
  tmpnames = names(list1)

  for(i in 1:length(list1)){
    v1 = list1[[tmpnames[i]]]
    matrix1 = list2[[tmpnames[i]]]
    commsample = intersect(names(v1),colnames(matrix1))
    v1 = v1[commsample]
    matrix1 = matrix1[,commsample]
    tmpcor = data.frame(cor = rep(0,nrow(matrix1)),
                        Pvalue = rep(1,nrow(matrix1)))
    rownames(tmpcor) = rownames(matrix1)
    for(k in 1:nrow(matrix1)){
      if(sum(!is.na(matrix1[k,])) < 3){next;}
      tmp = cor.test(v1, matrix1[k,])
      
      tmpcor$cor[k] = tmp$estimate
      tmpcor$Pvalue[k] = tmp$p.value
    }
    outcor[[i]] = tmpcor
  }
  names(outcor) = tmpnames
  return(outcor)
}

delete_dup_genes_forprotein <- function(expr,gtf.gene2symbol){
    # get most abundance genes when duplicated
    tnames = rownames(expr)
    gnames_org = gtf.gene2symbol[tnames,]$Gene
    sidx = sort.int(rowSums(expr,na.rm = T),decreasing = T,index.return = T)$ix
    tnames1 = tnames[sidx]
    gnames = gtf.gene2symbol[tnames1,]$Gene
    duplicateIDx = tnames1[duplicated(gnames)]
    idxx = !is.element(tnames,duplicateIDx)
    out = expr[idxx,]
    rownames(out) = gnames_org[idxx]
    out = out[rownames(out) != '',]
    return(out)
}

del_duplicate_rows_forpro <- function(tmpdata,pro.whole.nofilter.header){
  out = tmpdata
  out$gene = pro.whole.nofilter.header[rownames(out),]$Gene
  out = out[sort.int(tmpdata$Pvalue,decreasing = F,index.return = T)$ix,]
  out = out[!duplicated(out$gene) & !is.na(out$gene) & out$gene != '',]
  rownames(out) = out$gene
  return(out)
}


get_Metacorrlation_met2pro <- function(met.tissues,pro.tissues.v,Metapro,inputname){
  # Hypoxanthine
  teron = list()
  teron[[length(met.tissues)]] = ''
  #inputname = "Hypoxanthine"
  xnames = c()
  for(i in 1:length(met.tissues)){
    tmp = met.tissues[[i]]
    if(any(rownames(tmp)==inputname)){
      teron[[i]] = met.tissues[[i]][inputname,]
      xnames = c(xnames,names(met.tissues)[i])
    }
  }
  names(teron) = names(met.tissues)
  
  
  #get_correlation protein
  corTeron2pro = get_corr(teron,pro.tissues.v)
  corTeron2pro.Pvalue = str2num(list_element_select(corTeron2pro,rownames(Metapro),'Pvalue'))
  corTeron2pro.cor = str2num(list_element_select(corTeron2pro,rownames(Metapro),'cor'))
  corTeron2pro.Pvalue[is.na(corTeron2pro.Pvalue)] =1
  
  x<-list(p=corTeron2pro.Pvalue)
  xx.pro = MetaDE.pvalue(x = x,meta.method = 'AW') 
  tmpdata = data.frame(ID = rownames(corTeron2pro.Pvalue),
                       log2FC = rowMeans(corTeron2pro.cor,na.rm = T),
                       Pvalue = as.vector(xx.pro$meta.analysis$FDR))
  #tmpdata = del_duplicate_rows_forpro(tmpdata,pro.whole.nofilter.header)
  #tmpdata$ID = tmpdata$gene
  rownames(tmpdata) = tmpdata$ID
  tmpdata$Pvalue[tmpdata$Pvalue < 1e-30] = 1e-30
  return(tmpdata)
}

plot_heatmap_fromMetadata <- function(Metamrna,threecomm,filename = filename,tissue.systems = tissue.systems,
                                      height = 6,width = 8,cluster_cols = F,cluster_rows = T,
                                      show_rownames = T,
                                      show_colnames = T,
                                      fontsize_row = 11,fontsize_col = 11) {
  #mcc mrna
  mrnafc.common = Metamrna[threecomm, 
                           substr(colnames(Metamrna),1,7) == 'log2FC_']
  colnames(mrnafc.common) =gsub('log2FC_','',colnames(mrnafc.common))
  
  mrnapval.common = Metamrna[threecomm, 
                             substr(colnames(Metamrna),1,7) == 'Pvalue_']
  colnames(mrnapval.common) =gsub('Pvalue_','',colnames(mrnapval.common))
  
  
  mrnapval.common = mrnapval.common[rownames(mrnafc.common),colnames(mrnafc.common)]
  
  #
  display_matrix = matrix(' ',nrow(mrnapval.common),ncol(mrnapval.common))
  display_matrix[mrnapval.common < 0.1] = '.'
  display_matrix[mrnapval.common >= 0.01 & mrnapval.common < 0.05] = '*'
  display_matrix[mrnapval.common >= 0.001 & mrnapval.common < 0.01] = '**'
  display_matrix[mrnapval.common < 0.001] = '***'
  
  enbrks<-c(-1,-0.8,-0.4,-0.2,-0.05,-0.02,0.02,0.05,0.2,0.4,0.8,1)
  #DAscore.all_1 = DAscore.all
  #DAscore.all_1[DAscore.all_1 == 0] = NA
  
  cname = colnames(mrnafc.common)
  
  tclass = data.frame(System = tissue.systems[cname],row.names = cname)
  mrnafc.common.v = mrnafc.common
  mrnafc.common.v[mrnafc.common.v > 1] =1
  mrnafc.common.v[mrnafc.common.v < -1] = -1
  
  ann_colors = list(
    System = tissue.color)
  
  heatmap_mrna = pheatmap::pheatmap(mrnafc.common.v,scale = 'none',cluster_rows = cluster_rows,
                                    cluster_cols = cluster_cols, annotation_colors = ann_colors,
                                    fontsize_row = fontsize_row,fontsize_col = fontsize_col,
                                    annotation_col = tclass, show_rownames = show_rownames,
                                    show_colnames = show_rownames,
                                    breaks = enbrks,
                                    angle_col= 45,
                                    treeheight_row = 20,treeheight_col = 20,legend = T,
                                    display_numbers = display_matrix,
                                    #color=colorRampPalette(c('#3B4992','gray95','#BB0021'))(9),
                                    #colorRampPalette(c('#008280','gray95','#BB0021'))(11),
                                    color=colorRampPalette(c('#4DBBD5FF','gray95','#E64B35FF'))(11),
                                    #color=colorRampPalette(c('#3B4992','gray99','#EE0000'))(11),
                                    #color=colorRampPalette(c('#008280FF','gray99','#E64B35FF'))(11),
                                    file =filename,
                                    height = height,width = width)
  return(heatmap_mrna)
}

