library(genefilter)
library(limma)

# 1. General ##################

file2frame <- function(file, header = TRUE, sep = "\t",stringsAsFactors = F,row.names = NULL){
    out = read.delim(file, header = header, sep = sep,stringsAsFactors = stringsAsFactors,row.names = row.names)
    return(out)
}
writetxt <- function(x,filename,sep = '\t',row.names = F,col.names = T){
  write.table(x,filename,sep = sep,row.names = row.names,col.names = col.names,quote = F)
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

is.outliner <- function(x,coef = 3){
   out =  x %in% boxplot.stats(x,coef = coef)$out
    return(out)
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

summary_vector <- function(v){
  valid = v[!is.na(v)]
  ss = unique(valid)
  rr = rep(0,length(ss))
  tt = rep(0,length(ss))
  k = 1
  for (xx in ss){
    rr[k] = sum(valid == xx)
    k = k+1
  }
  tt = rr/length(valid) *100
  out = data.frame(Numbers = ceiling(rr),percentage = round(tt,2))
  rownames(out) = ss
  return(out)
}

pval_hyp <- function(n,gn,s,gs){
  out = 1-phyper(gs-1,gn,n-gn,s)
  return(out)
}


fillgaps_rowMedian <- function(X){
  out = X
  for(i in 1:nrow(out)){
    e = median(out[i,out[i,] > 1e-6])
    id = which(out[i,] < 1e-6)
    out[i,id] = e
  }
  return(out)
}

fillgaps_rowMin <- function(X){
  out = X
  for(i in 1:nrow(out)){
    e = min(out[i,out[i,] > 1e-6])
    id = which(out[i,] < 1e-6)
    out[i,id] = e
  }
  return(out)
}

list_to_matrix <- function(DEproFC){
    DEproFC_matrix = list()
    for(i in 1:length(DEproFC)){
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

#2. plot #####################
lghplot.boxplot <- function(p,width = 0.3){
  p = p +geom_violin(alpha = 0.5,outlier.size = 0) +
    geom_boxplot(width = width,color = 'black',fill = 'white',outlier.size = 0,alpha = 1)

}
lghplot.addtheme <- function (size = 18,hjust = FALSE,legend.position = "none"){
 p = theme(axis.text.x = element_text(size = size-2,face = 'bold',color = "black"))+
    theme(axis.text.y = element_text(size = size-2,face = 'bold',color = "black"))+
    theme(axis.title=element_text(size=size,face = 'bold',color = "black")) +
    theme(axis.title=element_text(size=size,face = 'bold',color = "black")) +
    theme(title = element_text(size= size,face = 'bold',color = "black"))+
    theme(#panel.grid.major =element_blank(),
        strip.text = element_text(size=size,face = 'bold',color = "black"),
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

lghplot.addthemeA <- function (size = 18,sizex = 18,sizey = 18,hjust = FALSE,legend.position = "none"){
 p = theme(axis.text.x = element_text(size = sizex,color = "black"))+
    theme(axis.text.y = element_text(size = sizey,color = "black"))+
    theme(axis.title=element_text(size=size,color = "black")) +
    theme(title = element_text(size= size,face = 'bold',color = "black"))+
    theme(legend.text=element_text(size=size-3,))+
    theme(#panel.grid.major =element_blank(),
        strip.text = element_text(size=size,color = "black"),
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

plot_DAscore <- function(DAscore,subsystem,fdr,
                         gtitle = 'DAscore Plot',
                         gcolors = c('grey','red','blue'),
                         gsize = 12,
                         xlabel = 'DAscore', gtype = 'fdr',
                         lengend.position = 'right'){
  idx = sort.int(DAscore,decreasing = F,index.return = T)$ix
  FDR = cut(-fdr,breaks = c(-1.1,-0.05,-0.01,-0.001,0.1),labels = c('N.S','<0.05','<0.01','<0.001'))
  
  
  # FDR = rep(1,length(DAscore))
  # FDR[fdr >= 0.01 & fdr < 0.05] = 2
  # FDR[fdr > 0.001 & fdr < 0.01] = 3
  # FDR[fdr <= 0.001] = 4
  # tlabel = c('N.S','<0.05','<0.01','<0.001')
  # FDR = factor(FDR,label = tlabel[unique(sort(FDR))])
  
  Significance = rep(1,length(DAscore))
  Significance[DAscore > 0 & fdr < 0.05] = 2
  Significance[DAscore < 0 & fdr < 0.05] = 3
  tlabel = c('N.S','UP','Down')
  tid = unique(sort(Significance))
  Significance = factor(Significance, label = tlabel[tid])
  
  plotdata = data.frame(DAscore = DAscore[idx],
                        subsystem = factor(subsystem[idx],level = subsystem[idx]),
                        fdr = fdr[idx],
                        Significance = Significance[idx],
                        FDR = FDR[idx])
  #ylab(bquote(~-Log[10]~italic(FDR)))
  if (gtype == 'fdr'){
      p = ggplot(plotdata,aes(x = subsystem,y = DAscore)) + scale_shape_discrete(solid=T)+
        geom_point(aes(colour = Significance,size = FDR))+
        scale_colour_manual(values=gcolors[tid])+ ggtitle(gtitle)+
        xlab('')+ ylab(xlabel)+ 
        coord_flip()+lghplot.addtheme(size =  gsize,legend.position = lengend.position)
  }else{
    plotdata$Pvalue = FDR[idx]
    p = ggplot(plotdata,aes(x = subsystem,y = DAscore)) + scale_shape_discrete(solid=T)+
      geom_point(aes(colour = Significance,size = Pvalue))+
      scale_colour_manual(values=gcolors[tid])+ ggtitle(gtitle)+
      xlab('')+ ylab(xlabel)+ 
      coord_flip()+lghplot.addtheme(size =  gsize,legend.position = lengend.position)
  }
  return(p)
}


# 3.  DEgenes/DEflux ##################

rowSd<-function (x, ...) {
    sqr = function(x) x * x
    n = rowSums(!is.na(x))
    n[n <= 1] = NA
    return(sqrt(rowSums(sqr(x - rowMeans(x, ...)), ...)/(n - 1)))
}

# calculate grpMedian
grpMeans <- function(expr, grp){
  grpList <- unique(grp)
  res <- data.frame(matrix(NA, nrow=nrow(expr), ncol=length(grpList)))
  colnames(res) <- grpList
  rownames(res) <- rownames(expr)
  for (i in 1:length(grpList)){
      res[,i]=rowMeans(expr[, grp==grpList[i]],na.rm = T)
  }
  return(res)
}

grpSds <- function(expr, grp){
  grpList <- unique(grp)
  res <- data.frame(matrix(NA, nrow=nrow(expr), ncol=length(grpList)))
  colnames(res) <- grpList
  rownames(res) <- rownames(expr)
  for (i in 1:length(grpList)){
      res[,i]=rowSds(expr[, grp==grpList[i]],na.rm = T)
  }
  return(res)
}

grpMeanSds <-function(expr, grp, ndigits=2){
  grpmeans=round(grpMeans(expr,grp), ndigits)
  grpsds=round(grpSds(expr, grp), ndigits)
  res <- sprintf("%s +/- %s", unlist(grpmeans), unlist(grpsds))
  res <- data.frame(matrix(res, nrow=nrow(grpmeans), byrow=F))
  colnames(res)=colnames(grpmeans)
  res$ID <- rownames(grpmeans)
  return(res)
}

DEGenes <- function(expr, design, contrast=NULL, bayes.prior=0.01, adjust.method="fdr", annot=NULL){
  require(limma)
  thenames = rownames(expr)
  if (is.null(contrast)){
      contrast=data.frame(matrix(0, nrow=ncol(design), ncol=ncol(design)))
      rownames(contrast)=colnames(contrast)=colnames(design)
  }
  fit1=lmFit(expr, design)
  fit1=contrasts.fit(fit1, contrast)
  fit1=eBayes(fit1, bayes.prior)
   
  tt1=topTable(fit1, coef=1, adjust.method=adjust.method, number=nrow(expr))
  tt1$ID=rownames(tt1)
  tt1=tt1[,c("ID", "AveExpr", "logFC", "P.Value", "adj.P.Val")]
  #fc.sign=rep(1, nrow(tt1))
  #fc.sign[sign(tt1$logFC)==-1]=-1
  #tt1$logFC=fc.sign*2^abs(tt1$logFC)
  colnames(tt1)[3:5]=sprintf("%s.%s", rep(colnames(contrast)[1],3), c("log2FC", "PV", "FDR"))
  res=tt1
  if (ncol(contrast)>1){
    for (i in 2:ncol(contrast)){
      tt1=topTable(fit1, coef=i, adjust.method="fdr", number=nrow(expr))
      tt1=tt1[,c("ID", "logFC", "P.Value", "adj.P.Val")]
      #fc.sign=rep(1, nrow(tt1))
      #fc.sign[sign(tt1$logFC)==-1]=-1
      #tt1$logFC=fc.sign*2^abs(tt1$logFC)
      colnames(tt1)[2:4]=sprintf("%s.%s", rep(colnames(contrast)[i],3), c("log2FC", "PV", "FDR"))
      res=merge(res, tt1, by="ID")
    }
  }
  if (!is.null(annot)) res=merge(annot, res, by="ID")
  grp=as.vector(apply(design, 1, function(x) paste(rep(colnames(design), x), sep="", collapse=".")))
  res=merge(res, grpMeanSds(expr, grp), by="ID")
  rownames(res) = res$ID
  res = res[thenames,]
  colnames(res) = c('ID', 'AveExpr', 'log2FC', 'Pvalue' ,'FDR', 'Case', 'Control')
  return(res)
}


DEGenes.simplified <-   function(m,catagory,subset = NULL){
  ## catagory: Case and Control
  if (is.logical(catagory)){
    cc = rep('Control',length(catagory))
    cc[catagory] = 'Case'
  }else{
    cc = catagory
  }
  if (!is.null(subset)){
     cc = cc[subset]
     m = m[,subset]
  }
  design = model.matrix(~-1+cc)
  colnames(design)=c("Case", "Control")
  contrast = makeContrasts(Case-Control, levels=colnames(design))
  res = DEGenes(m, design = design, contrast=contrast)
  colnames(res) = c('ID', 'AveExpr', 'log2FC', 'Pvalue' ,'FDR', 'Case', 'Control')
  return(res)
}


# 4. flux pathway analysis #################
Flux.subsystem_WDA <- function(DEflux,annote,pdffile = NULL,cutoff.nrxns = 1,cutoff.pval = 0.05,cutoff.fc = 0,omit.pathways = NULL){
  ax = summary(DEflux$log2FC)
  idxx = DEflux$Pvalue < cutoff.pval
  idup = DEflux$Pvalue < cutoff.pval & DEflux$log2FC > cutoff.fc #ax[[5]]
  iddown = DEflux$Pvalue < cutoff.pval & DEflux$log2FC < -cutoff.fc #ax[[2]]
  u_subsystem = unique(annote$subSystemes)
  DAout = data.frame(subsystem = u_subsystem,
                     DAscore = rep(0,length(u_subsystem)),
                     stringsAsFactors = F)
  annote_up = annote[idup,]
  annote_down = annote[iddown,]
  nsig = sum(idxx)
  nrxns = nrow(annote)
  for(i in 1:length(u_subsystem)){
    if (u_subsystem[i] == 'Transport, extracellular'){ next;}
    if (u_subsystem[i] == ''){ next;}
    if (sum(u_subsystem[i] == omit.pathways) > 0 ){next;}
    tmp_up = sum(annote_up$subSystemes == u_subsystem[i])
    tmp_down = sum(annote_down$subSystemes == u_subsystem[i])
    tmp_all = sum(annote$subSystemes == u_subsystem[i])
    if(tmp_all < cutoff.nrxns) { next;}
    #if(abs(tmp_up + tmp_down) < cutoff.nrxns) { next;}
    weight = 100*(tmp_up + tmp_down)/nsig
    DAout$DAscore[i] = weight*(tmp_up - tmp_down)/tmp_all
    #DAout$DAscore[i] = (tmp_up - tmp_down)/tmp_all
    #if(tmp_all < 2) {DAout$DAscore[i] =0 }
  }
  DAout = DAout[abs(DAout$DAscore) > 1e-3,]
  idx = sort.int(DAout$DAscore,decreasing = F,index.return =T )$ix
  DAout = DAout[idx,]
  
  bb=DAout$subsystem
  DAout$subsystem = factor(bb,level = bb)
  DAforplot  = data.frame(subsystem = factor(bb,level = bb),
                          WDAscore =  DAout$DAscore,
                          type=factor(c(rep("Down Pathway", sum(DAout$DAscore < 0)), 
                                        rep("Up Pathway", sum(DAout$DAscore > 0))), 
                                      level = c("Down Pathway","Up Pathway"))
  )
  
  CPCOLS <- c("#8DA1CB", "#FD8D62")
  p =  ggplot(data=DAforplot, aes(y=WDAscore, x=subsystem,fill = type)) +
    geom_bar(stat="identity", width=0.8) + theme(legend.position="none")+ theme_bw() + coord_flip()+
    scale_fill_manual(values = CPCOLS)+
    theme(axis.text=element_text(face = "bold", color="gray50")) +
    #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(axis.text.y = element_text(size = 10,face = 'bold'))
  print(p)
  if (!is.null(pdffile)){
    pdf(pdffile)
    print(p)
    dev.off()
  }
  return(DAforplot) 
}

Flux.subsystem_WDA_bootstrap <- function(DEflux,bootstrap_fc,bootstrap_pval,annote, cutoff.nrxns = 1,cutoff.pval = 0.05,cutoff.fc = 0,omit.pathways = NULL){
  ax = summary(DEflux$log2FC)
  idxx = DEflux$Pvalue < cutoff.pval
  idup = DEflux$Pvalue < cutoff.pval & DEflux$log2FC > cutoff.fc #ax[[5]]
  iddown = DEflux$Pvalue < cutoff.pval & DEflux$log2FC < -cutoff.fc #ax[[2]]
  u_subsystem = unique(annote$subSystemes)
  DAout = data.frame(subsystem = u_subsystem,
                     DAscore = rep(0,length(u_subsystem)),
		                 DApvalue = rep(1,length(u_subsystem)),
                     stringsAsFactors = F)
  annote_up = annote[idup,]
  annote_down = annote[iddown,]
  #bootstrap
  bootstrap_up = bootstrap_fc > cutoff.fc & bootstrap_pval < cutoff.pval
  bootstrap_down = bootstrap_fc < cutoff.fc & bootstrap_pval < cutoff.pval
  bootstrap_nsig = colSums(bootstrap_pval < cutoff.pval)

  nboots = ncol(bootstrap_fc)
  
  nsig = sum(idxx)
  nrxns = nrow(annote)
  for(i in 1:length(u_subsystem)){
    if (u_subsystem[i] == 'Transport, extracellular'){ next;}
    if (u_subsystem[i] == ''){ next;}
    if (sum(u_subsystem[i] == omit.pathways) > 0 ){next;}
    tmp_up = sum(annote_up$subSystemes == u_subsystem[i])
    tmp_down = sum(annote_down$subSystemes == u_subsystem[i])
    tmp_all = sum(annote$subSystemes == u_subsystem[i])
    if(tmp_all < cutoff.nrxns) { next;}
    #if(abs(tmp_up + tmp_down) < cutoff.nrxns) { next;}
    weight = 1#100*(tmp_up + tmp_down)/nsig
    #weight = 100*(tmp_up + tmp_down)/nsig
    DAout$DAscore[i] = weight*(tmp_up - tmp_down)/tmp_all
    #DAout$DAscore[i] = (tmp_up - tmp_down)/tmp_all
    #if(tmp_all < 2) {DAout$DAscore[i] =0 }
    
    # bootstrap
    idx = annote$subSystemes == u_subsystem[i]
    tmp_boots_up = colSums(bootstrap_up[idx,])
    tmp_boots_down = colSums(bootstrap_down[idx,])
    #boots_weight = 100*(tmp_boots_up + tmp_boots_down)/bootstrap_nsig
    boots_weight = 1#100*(tmp_boots_up + tmp_boots_down)/bootstrap_nsig
    boots_DA = weight*(tmp_boots_up - tmp_boots_down)/tmp_all
    if (DAout$DAscore[i] > 0){
      DAout$DApvalue[i] = sum(boots_DA > DAout$DAscore[i])/nboots
    }else if(DAout$DAscore[i] < 0){
      DAout$DApvalue[i] = sum(boots_DA < DAout$DAscore[i])/nboots
    }else{
      DAout$DApvalue[i] = 1
    }
  }
  DAout = DAout[abs(DAout$DAscore) > 1e-3,]
  idx = sort.int(DAout$DAscore,decreasing = F,index.return =T )$ix
  DAout = DAout[idx,]
  DAout$fdr = p.adjust(DAout$DApvalue,method = 'fdr')
  return(DAout) 
}






Flux.subsystem_WDA_bak <- function(DEflux,annote,pdffile = NULL,cutoff.nrxns = 1,cutoff.pval = 0.05,omit.pathways = NULL){
    idxx = DEflux$Pvalue < cutoff.pval
    idup = DEflux$Pvalue < cutoff.pval & DEflux$log2FC > 0
    iddown = DEflux$Pvalue < cutoff.pval & DEflux$log2FC < 0
    u_subsystem = unique(annote$subSystemes)
    DAout = data.frame(subsystem = u_subsystem,
                       DAscore = rep(0,length(u_subsystem)),
                       stringsAsFactors = F)
    annote_up = annote[idup,]
    annote_down = annote[iddown,]
    nsig = sum(idxx)
    nrxns = nrow(annote)
    for(i in 1:length(u_subsystem)){
        if (u_subsystem[i] == 'Transport, extracellular'){ next;}
        if (sum(u_subsystem[i] == omit.pathways) > 0 ){next;}
        tmp_up = sum(annote_up$subSystemes == u_subsystem[i])
        tmp_down = sum(annote_down$subSystemes == u_subsystem[i])
        tmp_all = sum(annote$subSystemes == u_subsystem[i])
	if(tmp_all < cutoff.nrxns) { next;}
        weight = 100*(tmp_up + tmp_down)/nsig
        DAout$DAscore[i] = weight*(tmp_up - tmp_down)/tmp_all
	#DAout$DAscore[i] = (tmp_up - tmp_down)/tmp_all
        #if(tmp_all < 2) {DAout$DAscore[i] =0 }
    }
    DAout = DAout[abs(DAout$DAscore) > 1e-6,]
    idx = sort.int(DAout$DAscore,decreasing = F,index.return =T )$ix
    DAout = DAout[idx,]

    bb=DAout$subsystem
    DAout$subsystem = factor(bb,level = bb)
    DAforplot  = data.frame(subsystem = factor(bb,level = bb),
                             WDAscore =  DAout$DAscore,
                             type=factor(c(rep("Down Pathway", sum(DAout$DAscore < 0)), 
                                           rep("Up Pathway", sum(DAout$DAscore > 0))), 
                                         level = c("Down Pathway","Up Pathway"))
                            )

     CPCOLS <- c("#8DA1CB", "#FD8D62")
     p =  ggplot(data=DAforplot, aes(y=WDAscore, x=subsystem,fill = type)) +
            geom_bar(stat="identity", width=0.8) + theme(legend.position="none")+ theme_bw() + coord_flip()+
            scale_fill_manual(values = CPCOLS)+
            theme(axis.text=element_text(face = "bold", color="gray50")) +
	    #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            theme(axis.text.y = element_text(size = 10,face = 'bold'))
     print(p)
     if (!is.null(pdffile)){
       pdf(pdffile)
       print(p)
       dev.off()
     }
     return(DAforplot) 
}

Flux.subsystem_enrichment <- function(DEflux,annote,direct){
  require(clusterProfiler)
  if(direct == 'UP'){
      idxx = DEflux$Pvalue < 0.05  & DEflux$log2FC > 0
   }else if(direct == 'DOWN'){
      idxx = DEflux$Pvalue < 0.05  & DEflux$log2FC < 0
   }else{
      idxx = DEflux$Pvalue < 0.05  #& DEflux$log2FC < 0
   }

  sigflux = rownames(DEflux)[idxx]
  fluxgmt = data.frame(ont = annote$subSystemes,gene = rownames(annote),stringsAsFactors = F) 
  
  Recon3D <- enricher(gene = sigflux,
                TERM2GENE=fluxgmt,
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =1,
                qvalueCutoff = 1
                )
  Recon3Dout = Recon3D@result
  Recon3Dout$DB = rep('Recon3D',dim(Recon3Dout)[1])
  n = dim(Recon3Dout)[2]
  Recon3Dout = Recon3Dout[,c(n,1:(n-1))]
  return(Recon3Dout)
}
#

#### flux colormap
showpanel <- function(col)
{
  image(z=matrix(1:100, ncol=1), col=col, xaxt="n", yaxt="n" )
}
hex2rgb <- function(x,sep = ", "){
  y = rep('',length(x))
  for(i in 1:length(x)){
     y[i] = paste(as.vector(col2rgb(x[i])), collapse = sep)
  }
  return(y)
}

flux.colormap <- function(insvgfile,outsvgfile,rxns,log2fc, col = NULL, brks = NULL){

  
  require(gplots)
  xx = readLines(insvgfile)
  Index.rxn = which(regexpr('class="reaction"',xx) > 0)
  end.rxn = which(regexpr('<g id="nodes">',xx) > 0)
  Index.rxn = c(Index.rxn, end.rxn)
  
   if (is.null(col)){
    col = colorpanel(7,"blue","gray80","red")
    maxvalue = max(abs(log2fc),na.rm = T)
    brks = 2*(1:(length(col)/2))*maxvalue/length(col)
    brks = c(-brks,-1e-6,1e-6,brks)
  }
  tcolor = as.vector(cut(log2fc,breaks = brks,labels = col))
  rgbs = hex2rgb(tcolor)
  
  for(i in 1:(length(Index.rxn)-1)){
    tblock = xx[Index.rxn[i]:Index.rxn[i+1]]
    #rxn = regmatches(tblock,regexpr('>[^<]+</text>',tblock))
    rxn = regmatches(tblock,regexpr('"visible">[^<]+</text>',tblock))
    rxn = substring(rxn,11,nchar(rxn) - 7)
    if (substr(rxn,1,2) == 'EX_'){
      rxn = paste0(substr(rxn,1,nchar(rxn)-2),'(e)')
    }
    tid = which(rxns == rxn)
    if(length(tid) ==0 ){next;}
    trgb = rgbs[tid]
    if (is.na(trgb)){next;}
    
    for(j in 1:length(tblock)){
      tblock[j] = gsub('rgb\\(255, 255, 255\\)','#NOT_CHANGE_THIS_COLOR',tblock[j])
      tblock[j] = gsub('rgb\\([^)]+\\)',paste0('rgb(',trgb,')'),tblock[j])
      tblock[j] = gsub('#NOT_CHANGE_THIS_COLOR','rgb\\(255, 255, 255\\)',tblock[j])
    }
    xx[Index.rxn[i]:Index.rxn[i+1]] = tblock
  }
  write(xx,file = outsvgfile)
}

