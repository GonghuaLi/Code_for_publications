loadRData <- function(fileName){
  #loads an RData file, and returns it
  bb = load(fileName)
  out = list()
  for(i in 1:length(bb)){
    out[[i]] = get(bb[i])
  }
  names(out) = bb
  return(out)
}

GOenrichment_symbol <- function(genes){
    geneids = IDconverter_local(genes, genetype = "SYMBOL")
    out = GOenrichment(unique(geneids$ENTREZID))
    return(out)
}

run_Mfuzz <- function(pro.whole.symbol,pro.whole.info,pro.c = 8){
    require(Mfuzz)
    pro.whole.mean = aggregate(t(pro.whole.symbol), by=list(pro.whole.info$stage), FUN=mean, na.rm = T)
    pro.whole.mean.t = t(pro.whole.mean)
    vid = rowSums(is.na(pro.whole.mean.t)) < 1
    pro.whole.mean.t = pro.whole.mean.t[vid,]
    pro.eset = new("ExpressionSet",exprs = pro.whole.mean.t)
    pro.eset.std = filter.std(pro.eset,min.std=0)
    pro.eset.std.stand  =  standardise(pro.eset.std)
    pro.m = 1.5#mestimate(pro.eset.std.stand)
    pro.class = mfuzz(pro.eset.std.stand, c = pro.c,m = pro.m)
    out = list()
    out$class = pro.class
    out$eset.std.stand = pro.eset.std.stand
    return(out)
}

standardise_matrix <- function (m) 
{
    out = m
    for (i in 1:nrow(m)) {
        out[i, ] <- (m[i, ] - mean(m[i, ], na.rm = TRUE))/sd(m[i,], na.rm = TRUE)
    }
    out = out[-1,]
    return(out)
}

standardise_1 <- function (eset) 
{
    data <- exprs(eset)
    for (i in 1:dim(data)[[1]]) {
        data[i, ] <- (data[i, ] - data[i, 1])/sd(data[i,], na.rm = TRUE)
    }
    exprs(eset) <- data
    eset
}

standardise_2 <- function (eset) 
{
    data <- exprs(eset)
    for (i in 1:dim(data)[[1]]) {
        data[i, ] <- (data[i, ] - data[i, 2])/sd(data[i,], na.rm = TRUE)
    }
    exprs(eset) <- data
    eset
}

run_Mfuzz_normal_by1 <- function(pro.whole.symbol,pro.whole.info,pro.c = 8){
    require(Mfuzz)
    vid = pro.whole.info$stage != '5'
    pro.whole.symbol = pro.whole.symbol[,vid]
    pro.whole.info = pro.whole.info[vid,]
    pro.whole.mean = aggregate(t(pro.whole.symbol), by=list(pro.whole.info$stage), FUN=mean, na.rm = T)
    pro.whole.mean.t = t(pro.whole.mean)
    vid = rowSums(is.na(pro.whole.mean.t)) < 1
    pro.whole.mean.t = pro.whole.mean.t[vid,]
    pro.eset = new("ExpressionSet",exprs = pro.whole.mean.t)
    pro.eset.std = filter.std(pro.eset,min.std=0)
    pro.eset.std.stand  =  standardise_1(pro.eset.std)
    pro.m = 1.5#mestimate(pro.eset.std.stand)
    #pro.m = mestimate(pro.eset.std.stand)
    pro.class = mfuzz(pro.eset.std.stand, c = pro.c,m = pro.m)
    out = list()
    out$class = pro.class
    out$eset.std.stand = pro.eset.std.stand
    out$info = pro.whole.info
    return(out)
}

run_Mfuzz_normal <- function(pro.whole.symbol,pro.whole.info,pro.c = 8){
    require(Mfuzz)
    vid = pro.whole.info$stage != '5'
    pro.whole.symbol = pro.whole.symbol[,vid]
    pro.whole.info = pro.whole.info[vid,]
    pro.whole.mean = aggregate(t(pro.whole.symbol), by=list(pro.whole.info$stage), FUN=mean, na.rm = T)
    pro.whole.mean.t = t(pro.whole.mean)
    vid = rowSums(is.na(pro.whole.mean.t)) < 1
    pro.whole.mean.t = pro.whole.mean.t[vid,]
    pro.eset = new("ExpressionSet",exprs = pro.whole.mean.t)
    pro.eset.std = filter.std(pro.eset,min.std=0)
    pro.eset.std.stand  =  standardise(pro.eset.std)
    pro.m = 1.5#mestimate(pro.eset.std.stand)
    #pro.m = mestimate(pro.eset.std.stand)
    pro.class = mfuzz(pro.eset.std.stand, c = pro.c,m = pro.m)
    out = list()
    #out$mfuzz = pro.class
    out$class = pro.class
    out$eset.std.stand = pro.eset.std.stand
    out$info = pro.whole.info
    return(out)
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

get_lmExpr <- function(expr.tissue,metadata.tissue){
  ngene = nrow(expr.tissue)
  lmExpr.age = data.frame(meanExpr = rowMeans(expr.tissue),
                          beta.age = rep(0,ngene),p.age = rep(1,ngene))
  for (i in 1:nrow(expr.tissue)){
    if(sum(!is.na(expr.tissue[i,metadata.tissue$type != 'Elderly_Treated'])) < 5){
        next;
    }
    #if (sum(expr.tissue[i,metadata.tissue$type != 'Elderly_Treated']  > 0) <3){
    #  next;
    #}
    bx = lm(expr.tissue[i,] ~ metadata.tissue$age,subset = metadata.tissue$type != 'Elderly_Treated')
    #bx = lm(met.expr[i,] ~ (clin.met$Class=='C') +  (clin.met$Class=='F1') + (clin.met$Gender == 'Female'))
    lmExpr.age$beta.age[i] = bx$coefficients[2]
    lmExpr.age$p.age[i] = car::Anova(bx)$`Pr(>F)`[1]  
  }
  rownames(lmExpr.age) = rownames(expr.tissue)
  #theader = headers[rownames(lmExpr.age),]
  #lmExpr.age = cbind(theader,lmExpr.age)
  return(lmExpr.age)
}

get_CPM_counts <- function(counts,log = T){
    cts <- counts
    y <- DGEList(cts)
    y <- calcNormFactors(y, method = "TMM")
    cpms <- cpm(y,log = log)
    return(cpms)
}

delete_dup_genes <- function(expr,gtf.gene2symbol){
    # get most abundance genes when duplicated
    tnames = rownames(expr)
    gnames_org = gtf.gene2symbol[tnames,]$gene_name
    sidx = sort.int(rowSums(expr),decreasing = T,index.return = T)$ix
    tnames1 = tnames[sidx]
    gnames = gtf.gene2symbol[tnames1,]$gene_name
    duplicateIDx = tnames1[duplicated(gnames)]
    idxx = !is.element(tnames,duplicateIDx)
    out = expr[idxx,]
    rownames(out) = gnames_org[idxx]
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

imputeMin <- function(x){
    y = x
    y[y< 0.1] = NA
    tt = apply(y,1,min,na.rm =T)
    for(i in 1:nrow(y)){
        y[i, is.na(y[i,])] = tt[i]
    }
    return(y)
}

imputeRF <- function(x){
    y = x
    y[y< 0.1] = NA
    y = randomForest::rfImpute(y)
    return(y)
}

imputeKNN <- function(x){
    y = x
    y[y< 0.1] = NA
    y = impute::impute.knn(y)
    return(y)
}

get_network <- function( xr,xp,isdiag = TRUE){
    if (isdiag){
        yr = melt(xr)
        yp = melt(xp)
        xx = t(apply(as.matrix(yr[,1:2]),1,sort))
        tnames = apply(xx,1,paste0,collapse = '___')
        idx = !duplicated(tnames) & xx[,1] != xx[,2] & abs(yr$value) > 0.6 & yp$value < 0.01
        out = yr$value[idx]
        yy = data.frame(rho = out,stringsAsFactors =F)
        rownames(yy) = tnames[idx]
    }else{
        yr = melt(xr)
        yp = melt(xp)
        tnames = apply(as.matrix(yr[,1:2]),1,paste0,collapse = '___')
        idx = abs(yr$value) > 0.6 & yp$value < 0.01
        out = yr$value[idx]
        yy = data.frame(rho = out,stringsAsFactors =F)
        rownames(yy) = tnames[idx]
    }
    return(as.data.frame(t(yy)))
   
}

get_DEexpr.treat <- function(expr.tissue,metadata.tissue,control = 4){
    ddtype = metadata.tissue$stage
    DEexpr.treat = DEGenes.simplified(expr.tissue,catagory = ddtype ==5, 
                                    subset = (ddtype ==5 | ddtype ==control))
  return(DEexpr.treat)
}



#enrich plot for figure 8
met.class_enrichment <- function(mets,annote){
  require(clusterProfiler)

  vmet = intersect(mets,rownames(annote))
  fluxgmt = data.frame(ont = annote$sub_class,gene = rownames(annote),stringsAsFactors = F) 
  
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