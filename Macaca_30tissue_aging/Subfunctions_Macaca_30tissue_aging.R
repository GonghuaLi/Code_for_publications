standardise_matrix <- function(m) {
  out <- m
  for (i in 1:nrow(m)) {
    out[i, ] <- (m[i, ] - mean(m[i, ], na.rm = TRUE)) / sd(m[i, ], na.rm = TRUE)
  }
  out <- out[-1, ]
  return(out)
}

standardise_matrix_1 <- function(m, center = 1) {
  out <- m
  for (i in 1:nrow(m)) {
    out[i, ] <- (m[i, ] - m[i, center]) / sd(m[i, ], na.rm = TRUE)
  }
  out <- out[-1, ]
  return(out)
}

rbind2 <- function(m1, m2) {
  m1f <- as.data.frame(m1)
  # rownames(m1f) = rownames(m1f)
  m2f <- as.data.frame(m2)
  # rownames(m2f) = rownames(m2f)
  out <- as.matrix(as.data.frame(rbindlist(list(m1f, m2f), fill = T)))
  rownames(out)  <- c(rownames(m1f), rownames(m2f))
  return(out)
}

file2frame <- function(file, header = TRUE, sep = "\t", stringsAsFactors = F,
                       row.names = NULL) {
  out <- read.delim(file, header = header, sep = sep, stringsAsFactors = stringsAsFactors,
    row.names = row.names)
  return(out)
}

writetxt <- function(x, filename, sep = "\t", row.names = F, col.names = T) {
  write.table(x, filename, sep = sep, row.names = row.names,
    col.names = col.names, quote = F)
}

writetxt_forGPMM <- function(x, filename) {
  name <- rownames(x)
  y <- cbind(as.data.frame(name), x)
  writetxt(y, filename)
}

lghplot.addtheme <- function(size = 18, sizex = 18, sizey = 18, hjust = FALSE, legend.position = "none") {
  p <- theme(axis.text.x = element_text(size = sizex, color = "black")) +
    theme(axis.text.y = element_text(size = sizey, color = "black")) +
    theme(axis.title = element_text(size = size, color = "black")) +
    theme(title = element_text(size = size, face = "bold",
      color = "black")) + theme(legend.text = element_text(size = size - 3, )) +
    theme(strip.text = element_text(size = size, color = "black"), legend.position = legend.position,
      panel.grid.minor = element_blank(), panel.background = element_blank(),
      axis.line = element_line(colour = "black"))
  if (hjust) {
    p <- p + theme(axis.text.x = element_text(angle = 45,
      hjust = 1))
  }
  return(p)
}

CVs <- function(x, margin = 1) {
  calc_CV <- function(a) {
    return(sd(a) / mean(a))
  }
  return(apply(x, margin, calc_CV))
}

summary_vector <- function(v) {
  valid <- v[!is.na(v)]
  ss <- unique(valid)
  rr <- rep(0, length(ss))
  tt <- rep(0, length(ss))
  k <- 1
  for (xx in ss) {
    rr[k] <- sum(valid == xx)
    k <- k + 1
  }
  tt <- rr / length(valid) * 100
  out <- data.frame(Numbers = ceiling(rr), percentage = round(tt,
    2))
  rownames(out) <- ss
  return(out)
}

is.outliner <- function(x, coef = 3) {
  out <- x %in% boxplot.stats(x, coef = coef)$out
  return(out)
}

data_summary <- function(data, varname, groupnames) {
  require(plyr)
  summary_func <- function(x, col) {
    c(mean = mean(x[[col]], na.rm = TRUE),
      sd = sd(x[[col]], na.rm = TRUE))
  }
  data_sum <- ddply(data, groupnames, .fun = summary_func,
    varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#
shorten_names_forGO1 <- function(x, n_word = 5, n_char = 30) {
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) >
    n_char)) {
    if (nchar(x) > n_char)
      x <- substr(x, 1, n_char)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x,
      " ")[[1]]), n_word)], collapse = " "), "...", sep = "")
    return(x)
  } else {
    return(x)
  }
}

parse_metascape <- function(tpath, gokegg = TRUE) {
  # tpath = './results/differTypes/differ_typeI_vs_typeII_pro_up/metascape_result.xlsx'
  my_data <- as.data.frame(readxl::read_excel(tpath, sheet = "Enrichment"))
  if (gokegg) {
    my_data <- my_data[my_data$Category == "GO Biological Processes" | my_data$Category == "KEGG Pathway", ]
  }
  idx <- regexpr("Summary", my_data$GroupID) < 0
  my_data <- my_data[idx, ]
  my_data <- my_data[!duplicated(my_data$GroupID), ]
  goup <- -my_data$`Log(q-value)`
  names(goup) <- my_data$Description
  goup <- sort(goup, decreasing = T)
  return(goup)
}
plot_enrich <- function(goup, godown, n_char = 60) {
  CPCOLS <- c("#BB0021", "#3B4992")
  aa <- sort(godown[1:min(10, length(godown))], decreasing = F)
  bb <- sort(goup[1:min(10, length(goup))], decreasing = F)
  tgoname <- capitalize(c(names(aa), names(bb)))
  tgoname[duplicated(tgoname)] <- paste0(" ", tgoname[duplicated(tgoname)])
  tgoname.short <- as.vector(sapply(tgoname, shorten_names_forGO1, n_char = n_char))
  tmp <- data.frame(log10pval = as.vector(c(aa, bb)),
    direction = factor(c(rep("Down", length(aa)), rep("Up", length(bb))),
      levels = c("Up", "Down")),
    goname =   factor(tgoname, levels = tgoname),
    showname = factor(tgoname.short, levels = tgoname.short)
  )
  # print(tmp)
  # tlab = (sapply(levels(tmp$goname), shorten_names_forGO))
  tlab <- tmp$showname
  tmp$number <- factor(1:nrow(tmp))
  p <- ggplot(data = tmp, aes(x = number, y = log10pval, fill = direction)) +
    geom_bar(stat = "identity", width = 0.8) + coord_flip() +
    scale_fill_manual(values = CPCOLS) + theme_classic() +
    scale_x_discrete(labels = tlab) + xlab("GO term") + theme(axis.text = element_text(face = "bold",
      color = "gray0")) + theme(axis.text.y = element_text(size = 8, face = "bold")) + lghplot.addtheme(size = 10) +
    ylab("-log10(p.adjust)") + xlab("")
  return(p)

}

fmt_dcimals <- function(decimals = 0) {
  function(x) format(x, nsmall = decimals, scientific = FALSE)
}

rowSds <- function(x, ...) {
  sqr <- function(x) x * x
  n <- rowSums(!is.na(x))
  n[n <= 1] <- NA
  return(sqrt(rowSums(sqr(x - rowMeans(x, ...)), ...) / (n -
    1)))
}

DEGenes <- function(expr, design, contrast = NULL, bayes.prior = 0.01, adjust.method = "fdr", annot = NULL) {
  require(limma)

  # calculate grpMedian
  grpMeans <- function(expr, grp) {
    grpList <- unique(grp)
    res <- data.frame(matrix(NA, nrow = nrow(expr), ncol = length(grpList)))
    colnames(res) <- grpList
    rownames(res) <- rownames(expr)
    for (i in 1:length(grpList)) {
      res[, i] <- rowMeans(expr[, grp == grpList[i]], na.rm = T)
    }
    return(res)
  }

  grpSds <- function(expr, grp) {
    grpList <- unique(grp)
    res <- data.frame(matrix(NA, nrow = nrow(expr), ncol = length(grpList)))
    colnames(res) <- grpList
    rownames(res) <- rownames(expr)
    for (i in 1:length(grpList)) {
      res[, i] <- rowSds(expr[, grp == grpList[i]], na.rm = T)
    }
    return(res)
  }

  grpMeanSds <- function(expr, grp, ndigits = 2) {
    grpmeans <- round(grpMeans(expr, grp), ndigits)
    grpsds <- round(grpSds(expr, grp), ndigits)
    res <- sprintf("%s +/- %s", unlist(grpmeans), unlist(grpsds))
    res <- data.frame(matrix(res, nrow = nrow(grpmeans), byrow = F))
    colnames(res) <- colnames(grpmeans)
    res$ID <- rownames(grpmeans)
    return(res)
  }

  thenames <- rownames(expr)
  if (is.null(contrast)) {
    contrast <- data.frame(matrix(0, nrow = ncol(design), ncol = ncol(design)))
    rownames(contrast) <- colnames(contrast) <- colnames(design)
  }
  fit1 <- lmFit(expr, design)
  fit1 <- contrasts.fit(fit1, contrast)
  fit1 <- eBayes(fit1, bayes.prior)

  tt1 <- topTable(fit1, coef = 1, adjust.method = adjust.method, number = nrow(expr))
  tt1$ID <- rownames(tt1)
  tt1 <- tt1[, c("ID", "AveExpr", "logFC", "P.Value", "adj.P.Val")]
  colnames(tt1)[3:5] <- sprintf("%s.%s", rep(colnames(contrast)[1], 3), c("log2FC", "PV", "FDR"))
  res <- tt1
  if (ncol(contrast) > 1) {
    for (i in 2:ncol(contrast)) {
      tt1 <- topTable(fit1, coef = i, adjust.method = "fdr", number = nrow(expr))
      tt1 <- tt1[, c("ID", "logFC", "P.Value", "adj.P.Val")]
      colnames(tt1)[2:4] <- sprintf("%s.%s", rep(colnames(contrast)[i], 3), c("log2FC", "PV", "FDR"))
      res <- merge(res, tt1, by = "ID")
    }
  }
  if (!is.null(annot)) res <- merge(annot, res, by = "ID")
  grp <- as.vector(apply(design, 1, function(x) paste(rep(colnames(design), x), sep = "", collapse = ".")))
  res <- merge(res, grpMeanSds(expr, grp), by = "ID")
  rownames(res) <- res$ID
  res <- res[thenames, ]
  colnames(res) <- c("ID", "AveExpr", "log2FC", "Pvalue", "FDR", "Case", "Control")
  return(res)
}


DEGenes.simplified <-   function(m, catagory, subset = NULL) {
  require(limma)
  ## catagory: Case and Control
  if (is.logical(catagory)) {
    cc <- rep("Control", length(catagory))
    cc[catagory] <- "Case"
  } else {
    cc <- catagory
  }
  if (!is.null(subset)) {
    cc <- cc[subset]
    m <- m[, subset]
  }
  design <- model.matrix(~ -1 + cc)
  colnames(design) <- c("Case", "Control")
  contrast <- makeContrasts(Case - Control, levels = colnames(design))
  res <- DEGenes(m, design = design, contrast = contrast)
  colnames(res) <- c("ID", "AveExpr", "log2FC", "Pvalue", "FDR", "Case", "Control")
  return(res)
}




plot_DEflux <- function(DEflux, x = "log2FC", y = "Pvalue",
                        pcutoff = 0.05, num.showlab = 10,
                        title = "DEflux", FCcutoff = 0, alpha = 0.3,
                        xlab = bquote(~ Log[2] ~ "fold change"), fixpointsize = NULL,
                        ylab = bquote(~ -Log[10] ~ italic(P)),
                        labSize = 3, legendLabSize = 5) {
  ## this function is from rGPMM package ligonghua@mail.kiz.ac.cn
  require(EnhancedVolcano)
  res2 <- DEflux
  res2$log2FC <- res2[[x]]
  res2$Pvalue <- res2[[y]]
  res2$log2FC[is.na(res2$log2FC)] <- 0
  res2$Pvalue[is.na(res2$Pvalue)] <- 1
  keyvals <- rep("gray50", nrow(res2))
  names(keyvals) <- rep("NS", nrow(res2))

  keyvals[which(res2$log2FC > FCcutoff & res2$Pvalue < pcutoff)] <- "Brown"
  names(keyvals)[which(res2$log2FC > FCcutoff & res2$Pvalue < pcutoff)] <- "UP"

  keyvals[which(res2$log2FC < -FCcutoff & res2$Pvalue < pcutoff)] <- "darkblue"
  names(keyvals)[which(res2$log2FC < -FCcutoff & res2$Pvalue < pcutoff)] <- "Down"


  if (!is.element("ID", colnames(res2))) {
    res2$ID <- rownames(res2)
  } else {
    rownames(res2) <-  res2$ID
  }
  thelab <- res2$ID

  uptmp <- res2[res2$log2FC > 0, ]
  downtmp <- res2[res2$log2FC < 0, ]

  xid <- sort.int(uptmp$Pvalue, decreasing = F, index.return = T)$ix
  uptmp <- uptmp[xid, ]
  if (nrow(uptmp) > num.showlab) {
    uplab <- uptmp$ID[1:num.showlab]
  } else {
    uplab <- uptmp$ID
  }

  xid <- sort.int(downtmp$Pvalue, decreasing = F, index.return = T)$ix
  downtmp <- downtmp[xid, ]
  if (nrow(uptmp) > num.showlab) {
    downlab <- downtmp$ID[1:num.showlab]
  } else {
    downlab <- downtmp$ID
  }

  ymax <- max(-log10(res2$Pvalue)) + 1

  thelab <- c(uplab, downlab)
  vindx <- res2[thelab, ]$Pvalue < pcutoff # is.element(res2$ID,thelab)
  thelab <- thelab[vindx]
  # thelab = thelab[res2$Pvalue[vindx] < pcutoff]
  # thelab[!is.element(thelab,c(uplab,downlab))] = NA
  all_label <- res2$ID
  all_label[!is.element(all_label, thelab)] <- NA
  if (is.null(fixpointsize)) {
    p <- EnhancedVolcano(res2,
      ylim = c(0, ymax),
      lab = all_label,
      # selectLab = thelab,
      # lab = thelab,
      x = "log2FC",
      y = "Pvalue",
      xlim = c(min(res2[["log2FC"]], na.rm = TRUE) - 0.1, max(res2[["log2FC"]], na.rm = TRUE) +
        0.1),
      # ylim = c(0, max(-log10(res2[['Pvalue']]), na.rm = TRUE) + 5),
      title = title,
      border = "full",
      titleLabSize = 12,
      FCcutoff = FCcutoff,
      cutoffLineWidth = 0,
      cutoffLineType = "blank",
      axisLabSize = 12,
      subtitle = NULL,
      cutoffLineCol = "white",
      gridlines.minor = F,
      gridlines.major = F,
      xlab = xlab,
      ylab = ylab,
      pCutoff = pcutoff,
      colCustom = keyvals,
      colAlpha = 4 / 5,
      # legendPosition = 'bottom',
      caption = NULL,
      legendPosition = "none",
      legendLabSize = 5,
      legendIconSize = 3,
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      pointSize = -alpha * log10(res2$Pvalue), labSize  = labSize,
      colConnectors = "black",
      max.overlaps = Inf)
  } else {
    p <- EnhancedVolcano(res2,
      ylim = c(0, ymax),
      lab = all_label,
      # selectLab = thelab,
      # lab = thelab,
      x = "log2FC",
      y = "Pvalue",
      xlim = c(min(res2[["log2FC"]], na.rm = TRUE) - 0.1, max(res2[["log2FC"]], na.rm = TRUE) +
        0.1),
      # ylim = c(0, max(-log10(res2[['Pvalue']]), na.rm = TRUE) + 5),
      title = title,
      border = "full",
      titleLabSize = 12,
      FCcutoff = FCcutoff,
      cutoffLineWidth = 0,
      cutoffLineType = "blank",
      axisLabSize = 12,
      caption = NULL,
      subtitle = NULL,
      cutoffLineCol = "white",
      gridlines.minor = F,
      gridlines.major = F,
      xlab = xlab,
      ylab = ylab,
      pCutoff = pcutoff,
      colCustom = keyvals,
      colAlpha = 4 / 5,
      # legendPosition = 'bottom',
      # direction = "y",
      legendPosition = "none",
      legendLabSize = 5,
      legendIconSize = 3,
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      pointSize = fixpointsize, labSize  = labSize,
      colConnectors = "black",
      max.overlaps = Inf)
  }


  return(p)
}




list_to_matrix <- function(DEproFC, alltissues) {
  DEproFC_matrix <- list()
  for (i in 1:length(alltissues)) {
    tmp <- matrix(DEproFC[[i]], 1, length(DEproFC[[i]]))
    tmp <- as.data.frame(tmp)
    colnames(tmp) <- names(DEproFC[[i]])
    DEproFC_matrix[[i]] <- tmp
  }
  DEproFC_matrix <- t(as.matrix(rbindlist(DEproFC_matrix, fill = T)))
  colnames(DEproFC_matrix) <- names(DEproFC)
  # vid = rowSums(is.na(DEproFC_matrix)) < ncol(DEproFC_matrix)/2
  # DEproFC_matrix = DEproFC_matrix[vid,]
  return(DEproFC_matrix)
}

GOenrichment <- function(vgenes, OrgDb = "org.Hs.eg.db", organism = "hsa",
                         pvalueCutoff = 0.01, qvalueCutoff = 0.1, seed = 2025520) {
  set.seed(seed)
  get_ratio <- function(gene_ratio, bg_ratio) {
    as.numeric(limma::strsplit2(gene_ratio, split = "/"))
    gene_value <- as.numeric(limma::strsplit2(gene_ratio, "/")[, 1]) / as.numeric(limma::strsplit2(gene_ratio, "/")[, 2])
    bg_value <- as.numeric(limma::strsplit2(bg_ratio, "/")[, 1]) / as.numeric(limma::strsplit2(bg_ratio, "/")[, 2])
    return (gene_value / bg_value)
  }

  gene.ent <- clusterProfiler::bitr(vgenes, fromType = "SYMBOL",
    toType = c("ENTREZID"), OrgDb = OrgDb)
  target.gene <- unlist(gene.ent[, c("ENTREZID")])

  kegg <- clusterProfiler::enrichKEGG(gene = target.gene, organism = organism,
    pAdjustMethod = "BH",
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff)
  kegg <- clusterProfiler::setReadable(kegg, OrgDb = OrgDb, keyType = "ENTREZID")@result
  kegg$ratio <- get_ratio(kegg$GeneRatio, kegg$BgRatio)
  kegg$type <- rep("KEGG", nrow(kegg))

  bp <- clusterProfiler::enrichGO(gene = target.gene,
    keyType = c("ENTREZID"), OrgDb = OrgDb, ont = "BP",
    pAdjustMethod = "BH", pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff,
    readable = TRUE)
  bp <- clusterProfiler::setReadable(bp, OrgDb = OrgDb, keyType = "ENTREZID")@result
  bp$ratio <- get_ratio(bp$GeneRatio, bp$BgRatio)
  bp$type <- rep("BP", nrow(bp))

  out <- rbind(kegg, bp)
  return(out)
}

GOenrichment.C2 <- function(vgenes, TERM2GENE, seed = 2025520) {
  set.seed(seed)
  get_ratio <- function(gene_ratio, bg_ratio) {
    as.numeric(limma::strsplit2(gene_ratio, split = "/"))
    gene_value <- as.numeric(limma::strsplit2(gene_ratio, "/")[, 1]) / as.numeric(limma::strsplit2(gene_ratio, "/")[, 2])
    bg_value <- as.numeric(limma::strsplit2(bg_ratio, "/")[, 1]) / as.numeric(limma::strsplit2(bg_ratio, "/")[, 2])
    return (gene_value / bg_value)
  }

  em <- enricher(vgenes, TERM2GENE = m_t2g)@result
  em$ratio <- get_ratio(em$GeneRatio, em$BgRatio)
  em$type <- rep("C2", nrow(em))
  return(em)
}



met.class_enrichment <- function(mets, annote) {
  require(clusterProfiler)

  vmet <- intersect(mets, rownames(annote))
  fluxgmt <- data.frame(ont = annote$sub_class,
    gene = rownames(annote), stringsAsFactors = F)

  Recon3D <- enricher(gene = vmet,
    TERM2GENE = fluxgmt,
    pAdjustMethod = "BH",
    minGSSize = 1,
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  Recon3Dout <- Recon3D@result
  Recon3Dout$DB <- rep("HMDBclass", dim(Recon3Dout)[1])
  n <- dim(Recon3Dout)[2]
  Recon3Dout <- Recon3Dout[, c(n, 1:(n - 1))]
  return(Recon3Dout)
}

fillgaps_rowMinNA_1 <- function(X) {
  out <- X
  min_e <- min(out, na.rm = T)
  alpha <- max(min_e, 1)
  for (i in 1:nrow(out)) {
    id <- which(is.na(out[i, ]))
    if (!any(id)) {
      next
    }
    e <- min(out[i, ], na.rm = T)
    if (is.infinite(e)) {
      e <- min_e
    }
    # out[i,id] = e  + alpha*runif(n=length(id), min = -0.01, max=0.01)
    out[i, id] <- e
  }
  out[is.infinite(out)] <- NA
  return(out)
}

fillgaps_MinNA <- function(X) {
  out <- X
  min_e <- min(out, na.rm = T)
  alpha <- max(min_e, 1)
  for (i in 1:nrow(out)) {
    id <- which(is.na(out[i, ]))
    if (!any(id)) {
      next
    }
    # e = min(out[i,],na.rm = T)
    # if(is.infinite(e)){
    out[i, id] <- min_e * (1 + runif(n = length(id), min = -0.01, max = 0.01))
    # }else{
    #  out[i,id] = e
    # }

  }
  return(out)
}

fillgaps_MinNA_formeta <- function(X, nas, nascutoff = 0.3) {
  out <- X
  min_e <- min(out, na.rm = T)
  # alpha = max(min_e,1)
  for (i in 1:nrow(out)) {
    id <- which(is.na(out[i, ]))
    if (!any(id)) {
      next
    }
    # e = min(out[i,],na.rm = T)
    if (nas[i] < nascutoff) {
      # e = min(out[i,],na.rm = T)
      out[i, id] <- min_e * (1 + runif(n = length(id), min = -0.01, max = 0.01))
    } else {
      out[i, ] <- NA
    }
  }
  return(out)
}


fillgaps_MinNA_formeta_v1 <- function(X, nas, nascutoff = 0.3) {
  out <- X
  min_e <- min(out, na.rm = T)
  # alpha = max(min_e,1)
  for (i in 1:nrow(out)) {
    id <- which(is.na(out[i, ]))
    if (!any(id)) {
      next
    }
    # e = min(out[i,],na.rm = T)
    e <- min(out[i, ], na.rm = T)
    if (is.infinite(e)) {
      e <- min_e
    }
    if (nas[i] < nascutoff) {
      # e = min(out[i,],na.rm = T)
      # out[i,id] = min_e*(1+runif(n=length(id), min = -0.01, max=0.01))
      out[i, id] <- e
    } else {
      out[i, ] <- NA
    }
  }
  return(out)
}

str2num <- function(x) {
  rname <- rownames(x)
  cname <- colnames(x)
  y <- matrix(as.numeric(x), ncol = ncol(x))
  rownames(y) <- rname
  colnames(y) <- cname
  return(y)
}

loadRData <- function(fileName) {
  bb <- load(fileName)
  out_for_this_output_unique <- list()
  for (i in 1:length(bb)) {
    out_for_this_output_unique[[i]] <- get(bb[i])
  }
  names(out_for_this_output_unique) <- bb
  return(out_for_this_output_unique)
}


list_element_select <- function(x, rname, cname) {
  out <- matrix("", length(rname), length(x))
  for (i in 1:length(x)) {
    out[, i] <- unlist(x[[i]][rname, cname])
  }

  colnames(out) <- names(x)
  rownames(out) <- rname
  return(out)
}


get_corr <- function(list1, list2) {
  outcor <- list()

  tmpnames <- names(list1)

  for (i in 1:length(list1)) {
    v1 <- list1[[tmpnames[i]]]
    matrix1 <- list2[[tmpnames[i]]]
    commsample <- intersect(names(v1), colnames(matrix1))
    v1 <- v1[commsample]
    matrix1 <- matrix1[, commsample]
    tmpcor <- data.frame(cor = rep(0, nrow(matrix1)),
      Pvalue = rep(1, nrow(matrix1)))
    rownames(tmpcor) <- rownames(matrix1)
    for (k in 1:nrow(matrix1)) {
      if (sum(!is.na(matrix1[k, ])) < 3) {
        next
      }
      tmp <- cor.test(v1, matrix1[k, ])

      tmpcor$cor[k] <- tmp$estimate
      tmpcor$Pvalue[k] <- tmp$p.value
    }
    outcor[[i]] <- tmpcor
  }
  names(outcor) <- tmpnames
  return(outcor)
}

delete_dup_genes_forprotein <- function(expr, gtf.gene2symbol) {
  # get most abundance genes when duplicated
  tnames <- rownames(expr)
  gnames_org <- gtf.gene2symbol[tnames, ]$Gene
  sidx <- sort.int(rowSums(expr, na.rm = T), decreasing = T, index.return = T)$ix
  tnames1 <- tnames[sidx]
  gnames <- gtf.gene2symbol[tnames1, ]$Gene
  duplicateIDx <- tnames1[duplicated(gnames)]
  idxx <- !is.element(tnames, duplicateIDx)
  out <- expr[idxx, ]
  rownames(out) <- gnames_org[idxx]
  out <- out[rownames(out) != "", ]
  return(out)
}

del_duplicate_rows_forpro <- function(tmpdata, pro.whole.nofilter.header) {
  out <- tmpdata
  out$gene <- pro.whole.nofilter.header[rownames(out), ]$Gene
  out <- out[sort.int(tmpdata$Pvalue, decreasing = F, index.return = T)$ix, ]
  out <- out[!duplicated(out$gene) & !is.na(out$gene) & out$gene != "", ]
  rownames(out) <- out$gene
  return(out)
}

#
get_CPM_counts <- function(cts, log = T) {
  require(edgeR)
  # cts <- txi$counts
  y <- DGEList(cts)
  y <- calcNormFactors(y, method = "TMM")
  cpms <- cpm(y, log = log)
  return(cpms)
}

get_CPM_nonormal <- function(x) {
  total_reads_per_sample <- colSums(x)
  cpm_matrix <- sweep(x, 2, total_reads_per_sample, "/") * 1e6
  return(cpm_matrix)
}

delete_dup_genes <- function(expr, gtf.gene2symbol) {
  # get most abundance genes when duplicated
  tnames <- rownames(expr)
  gnames_org <- gtf.gene2symbol[tnames, ]$gene_name
  sidx <- sort.int(rowSums(expr), decreasing = T, index.return = T)$ix
  tnames1 <- tnames[sidx]
  gnames <- gtf.gene2symbol[tnames1, ]$gene_name
  duplicateIDx <- tnames1[duplicated(gnames)]
  idxx <- !is.element(tnames, duplicateIDx)
  out <- expr[idxx, ]
  rownames(out) <- gnames_org[idxx]
  return(out)
}

updateRNA <- function(rnaData, headers, tissue.systems, macacageneInfo, outpath) {
  # filter with Rin > 6 and mapped counts > 10M
  idx <- colSums(rnaData$rawCounts.whole) >= 10000000 & rnaData$rawCounts.whole.info$RIN >= 6
  rnaData$rawCounts.whole <- rnaData$rawCounts.whole[, idx]
  rnaData$rawCounts.whole.info <- rnaData$rawCounts.whole.info[idx, ]

  outdata <- rnaData

  tcounts <- delete_dup_genes(rnaData$rawCounts.whole, headers)
  tcounts <- tcounts[rownames(tcounts) != "", ]

  # consider coding genes
  tmptype <- macacageneInfo[rownames(tcounts), ]$Type
  idx <- tmptype == "protein_coding" | substr(tmptype, nchar(tmptype) - 4, nchar(tmptype)) == "_gene"
  idx[is.na(idx)] <- FALSE
  tcounts <- tcounts[idx, ]
  tcpm <- get_CPM_counts(tcounts)
  tcpm.nonormal <- get_CPM_nonormal(tcounts)

  batch <- rnaData$rawCounts.whole.info$batch
  RINs <- as.numeric(rnaData$rawCounts.whole.info$RIN)

  #
  metadata <- data.frame(batch = batch,
    RIN = RINs,
    log10counts = log10(colSums(rnaData$rawCounts.whole)),
    age = rnaData$rawCounts.whole.info$age,
    tissue = rnaData$rawCounts.whole.info$tissues)
  metadata$batch <- factor(metadata$batch)
  metadata$tissue <- factor(metadata$tissue)

  # Apply linear model to remove RINes
  # tcpm.batchad_pre = limma::removeBatchEffect(tcpm,covariates = metadata$RIN,
  #                                            design=model.matrix(~ age+tissue, data = metadata))
  # tcpm.batchad_pre = limma::removeBatchEffect(tcpm,covariates = as.matrix(metadata[,c("RIN","log10counts")]))
  tcpm.batchad_model2 <- limma::removeBatchEffect(tcpm, covariates = RINs, batch = batch)
  tcpm.batchad_pre <- limma::removeBatchEffect(tcpm, covariates = RINs)

  # Apply ComBat to remove batch effects
  tcpm.batchadj <- sva::ComBat(
    dat = tcpm.batchad_pre,
    mod = model.matrix(~ age + tissue, data = metadata),
    batch = batch,
    par.prior = TRUE,    # Use parametric adjustments
    prior.plots = FALSE  # Disable plotting of priors
  )

  # for tsne plot Figure 1
  mrna.whole.std <- standardise_matrix(tcpm)
  mrna.whole.std.adj <- standardise_matrix(tcpm.batchadj)
  # mrna.whole.std.adj1 = standardise_matrix(tcpm.batchadj1)

  ## PCA before and after
  # because the expression data  is locate in the same magnitude,
  pca_before <- prcomp(t(tcpm), cor = F)
  percPCA_before <- 100 * summary(pca_before)$importance

  pca_after <- prcomp(t(tcpm.batchadj), cor = F)
  percPCA_after <- 100 * summary(pca_after)$importance

  ## boxplot
  png(paste0(outpath, "/boxplot_mRNA_before_adj.png"), width = 2000, height = 500)
  boxplot(tcpm, ylab = "Normalized log2 CPM", cex = 0.2,
    xlab = paste0("Samples of all mRNAs", " (n= ", ncol(tcpm), ")")
  )
  dev.off()

  ## boxplot
  png(paste0(outpath, "/boxplot_mRNA_after_adj.png"), width = 2000, height = 500)
  boxplot(tcpm.batchadj, ylab = "Normalized log2 CPM", cex = 0.2,
    xlab = paste0("Samples of all mRNAs after batch adj", " (n= ", ncol(tcpm.batchadj), ")")
  )
  dev.off()


  # plot PC1 PC2 and batch
  pdf(paste0(outpath, "/CPM_PC1_vs_PC2_before_batch_adj.pdf"), width = 9, height = 7)
  pp <- ggplot(, aes(pca_before$x[, 1], pca_before$x[, 2], color = batch)) + geom_point() +
    theme_bw() + lghplot.addtheme(legend.position = "right") +
    xlab(paste("PC1(", as.character(round(percPCA_before[2, 1], 1)), "%)", sep = "")) +
    ylab(paste("PC2(", as.character(round(percPCA_before[2, 2], 1)), "%)", sep = "")) +
    labs(title = paste0("PC1 vs PC2 before batch adjusted in mRNA"))
  print(pp)
  dev.off()

  # plot PC1 PC2  after batch
  pdf(paste0(outpath, "/CPM_PC1_vs_PC2_after_batch_adj.pdf"), width = 9, height = 7)
  pp <- ggplot(, aes(pca_after$x[, 1], pca_after$x[, 2], color = batch)) + geom_point() +
    theme_bw() + lghplot.addtheme(legend.position = "right") +
    xlab(paste("PC1(", as.character(round(percPCA_after[2, 1], 1)), "%)", sep = "")) +
    ylab(paste("PC2(", as.character(round(percPCA_after[2, 2], 1)), "%)", sep = "")) +
    labs(title = paste0("PC1 vs PC2 after batch adjusted in mRNA"))
  print(pp)
  dev.off()


  ## plot counts vs PC1 or PC2
  pdf(paste0(outpath, "/PC1_vs_counts_before_adj.pdf"), width = 9, height = 7)
  tmpxx <- cor.test(pca_before$x[, 1], metadata$log10counts)
  tmptext <- paste0("r=", signif(tmpxx$estimate, 2), "; p=", signif(tmpxx$p.value, 2))
  pp <- ggplot(, aes(pca_before$x[, 1], metadata$log10counts, color = batch)) +
    geom_point() + scale_color_npg() +
    theme_bw() + lghplot.addtheme(legend.position = "right") + # stat_ellipse(lwd=1,level = 0.99) +
    annotate(geom = "text", x = summary(pca_before$x[, 1])[5],
      y = summary(RINs)[5], label = tmptext,
      color = "black", size = 6) +
    xlab(paste("PC1(", as.character(round(percPCA_before[2, 1], 1)), "%)", sep = "")) +
    ylab("Number of total counts(log10)")
  print(pp)
  dev.off()

  pdf(paste0(outpath, "/PC2_vs_counts_before_adj.pdf"), width = 9, height = 7)
  tmpxx <- cor.test(pca_before$x[, 2], metadata$log10counts)
  tmptext <- paste0("r=", signif(tmpxx$estimate, 2), "; p=", signif(tmpxx$p.value, 2))
  pp <- ggplot(, aes(pca_before$x[, 2], metadata$log10counts, color = batch)) +
    geom_point() + scale_color_npg() +
    theme_bw() + lghplot.addtheme(legend.position = "right") +
    annotate(geom = "text", x = summary(pca_before$x[, 2])[5],
      y = summary(RINs)[5], label = tmptext,
      color = "black", size = 6) +
    xlab(paste("PC2(", as.character(round(percPCA_before[2, 2], 1)), "%)", sep = "")) +
    ylab("Number of total counts(log10)")
  print(pp)
  dev.off()


  ## plot RINs vs PC1 or PC2
  pdf(paste0(outpath, "/PC1_vs_RIN_before_adj.pdf"), width = 9, height = 7)
  tmpxx <- cor.test(pca_before$x[, 1], RINs)
  tmptext <- paste0("r=", signif(tmpxx$estimate, 2), "; p=", signif(tmpxx$p.value, 2))
  pp <- ggplot(, aes(pca_before$x[, 1], RINs, color = batch)) +
    geom_point() + scale_color_npg() +
    theme_bw() + lghplot.addtheme(legend.position = "right") + # stat_ellipse(lwd=1,level = 0.99) +
    annotate(geom = "text", x = summary(pca_before$x[, 1])[5],
      y = summary(RINs)[5], label = tmptext,
      color = "black", size = 6) +
    xlab(paste("PC1(", as.character(round(percPCA_before[2, 1], 1)), "%)", sep = "")) +
    ylab("RIN value")
  print(pp)
  dev.off()

  pdf(paste0(outpath, "/PC2_vs_RIN_before_adj.pdf"), width = 9, height = 7)
  tmpxx <- cor.test(pca_before$x[, 2], RINs)
  tmptext <- paste0("r=", signif(tmpxx$estimate, 2), "; p=", signif(tmpxx$p.value, 2))
  pp <- ggplot(, aes(pca_before$x[, 2], RINs, color = batch)) +
    geom_point() + scale_color_npg() +
    theme_bw() + lghplot.addtheme(legend.position = "right") +
    annotate(geom = "text", x = summary(pca_before$x[, 2])[5],
      y = summary(RINs)[5], label = tmptext,
      color = "black", size = 6) +
    xlab(paste("PC2(", as.character(round(percPCA_before[2, 2], 1)), "%)", sep = "")) +
    ylab("RIN value")
  print(pp)
  dev.off()

  ## plot RIN vs PC1 or PC2 after batch adj
  pdf(paste0(outpath, "/PC1_vs_RIN_after_adj.pdf"), width = 9, height = 7)
  tmpxx <- cor.test(pca_after$x[, 1], RINs)
  tmptext <- paste0("r=", signif(tmpxx$estimate, 2), "; p=", signif(tmpxx$p.value, 2))
  pp <- ggplot(, aes(pca_after$x[, 1], RINs, color = batch)) +
    geom_point() + scale_color_npg() +
    theme_bw() + lghplot.addtheme(legend.position = "right") +
    annotate(geom = "text", x = summary(pca_after$x[, 1])[5],
      y = summary(RINs)[5], label = tmptext,
      color = "black", size = 6) +
    xlab(paste("PC1(", as.character(round(percPCA_after[2, 1], 1)), "%)", sep = "")) +
    ylab("RIN value")
  print(pp)
  dev.off()

  pdf(paste0(outpath, "/PC2_vs_RIN_after_adj.pdf"), width = 9, height = 7)
  tmpxx <- cor.test(pca_after$x[, 2], RINs)
  tmptext <- paste0("r=", signif(tmpxx$estimate, 2), "; p=", signif(tmpxx$p.value, 2))
  pp <- ggplot(, aes(pca_after$x[, 2], RINs, color = batch)) +
    geom_point() + scale_color_npg() +
    theme_bw() + lghplot.addtheme(legend.position = "right") +
    annotate(geom = "text", x = summary(pca_after$x[, 2])[5],
      y = summary(RINs)[5], label = tmptext,
      color = "black", size = 6) +
    xlab(paste("PC2(", as.character(round(percPCA_after[2, 2], 1)), "%)", sep = "")) +
    ylab("RIN value")
  print(pp)
  dev.off()


  ## plot counts vs PC1 or PC2 after batch adj
  pdf(paste0(outpath, "/PC1_vs_counts_after_adj.pdf"), width = 9, height = 7)
  tmpxx <- cor.test(pca_after$x[, 1], metadata$log10counts)
  tmptext <- paste0("r=", signif(tmpxx$estimate, 2), "; p=", signif(tmpxx$p.value, 2))
  pp <- ggplot(, aes(pca_after$x[, 1], metadata$log10counts, color = batch)) +
    geom_point() + scale_color_npg() +
    theme_bw() + lghplot.addtheme(legend.position = "right") + # stat_ellipse(lwd=1,level = 0.99) +
    annotate(geom = "text", x = summary(pca_after$x[, 1])[5],
      y = summary(RINs)[5], label = tmptext,
      color = "black", size = 6) +
    xlab(paste("PC1(", as.character(round(percPCA_before[2, 1], 1)), "%)", sep = "")) +
    ylab("Number of total counts(log10)")
  print(pp)
  dev.off()

  pdf(paste0(outpath, "/PC2_vs_counts_after_adj.pdf"), width = 9, height = 7)
  tmpxx <- cor.test(pca_after$x[, 2], metadata$log10counts)
  tmptext <- paste0("r=", signif(tmpxx$estimate, 2), "; p=", signif(tmpxx$p.value, 2))
  pp <- ggplot(, aes(pca_after$x[, 2], metadata$log10counts, color = batch)) +
    geom_point() + scale_color_npg() +
    theme_bw() + lghplot.addtheme(legend.position = "right") +
    annotate(geom = "text", x = summary(pca_after$x[, 2])[5],
      y = summary(RINs)[5], label = tmptext,
      color = "black", size = 6) +
    xlab(paste("PC2(", as.character(round(percPCA_before[2, 2], 1)), "%)", sep = "")) +
    ylab("Number of total counts(log10)")
  print(pp)
  dev.off()



  #
  ### for tsne for batch before adjust
  p1 <- tsne(mrna.whole.std, labels = tissue.systems[rnaData$rawCounts.whole.info$tissues],
    legendtextsize = 10, dotsize = 2, seed = 2025520)
  p1 <- p1 + theme_classic() + lghplot.addtheme(legend.position = "none") +
    scale_color_npg() + xlab("tsne 1") + ylab("tsne 2") +
    theme(axis.line = element_line(size = 1.0)) + ggtitle("Transcriptome")
  pdf(file = paste0(outpath, "/tSNE_mrna_using_standardised_before_adjust.pdf"), height = 4, width = 4)
  print(p1)
  dev.off()


  ### for tsne for batch after ajust
  p2 <- tsne(mrna.whole.std.adj, labels = tissue.systems[rnaData$rawCounts.whole.info$tissues],
    legendtextsize = 10, dotsize = 2, seed = 2025520)
  p2 <- p2 + theme_classic() + lghplot.addtheme(legend.position = "none") +
    scale_color_npg() + xlab("tsne 1") + ylab("tsne 2") +
    theme(axis.line = element_line(size = 1.0)) + ggtitle("Transcriptome")
  pdf(file = paste0(outpath, "/tSNE_mrna_using_standardised_after_adjust.pdf"), height = 4, width = 4)
  print(p2)
  dev.off()

  outdata$mrna.whole <- tcpm.batchadj
  outdata$mrna.whole.info <- outdata$rawCounts.whole.info

  # tissues
  tmp.mrna.tissues <- list()
  tmp.mrna.tissues.info <- list()
  tmp.mrna.tissues.org <- list()
  tmp.mrna.tissues.adjRIN <- list()
  tmp.mrna.tissues.org.info <- list()
  tmp.mrna.tissues.model2 <- list()
  tmp.tissue.names <- unique(outdata$rawCounts.whole.info$tissues)
  for (i in 1:length(tmp.tissue.names)) {
    thisname <- tmp.tissue.names[i]
    vids <- outdata$rawCounts.whole.info$tissues == thisname
    # get high confidance tissue data
    thiscouts <- tcounts[, vids]
    tmpcpm.nonormal <- tcpm.nonormal[, vids]
    xids <- rowMeans(tmpcpm.nonormal) > 1 & rowSums(thiscouts < 5) < ncol(thiscouts) * 0.2
    tmp.mrna.tissues[[thisname]] <- tcpm.batchadj[xids, vids]
    tmp.mrna.tissues.info[[thisname]] <- outdata$rawCounts.whole.info[vids, ]
    tmp.mrna.tissues.adjRIN[[thisname]] <- tcpm.batchad_pre[xids, vids]
    tmp.mrna.tissues.model2[[thisname]] <- tcpm.batchad_model2[xids, vids]
    tmp.mrna.tissues.org[[thisname]] <- tcpm[xids, vids]
    tmp.mrna.tissues.org.info[[thisname]] <- outdata$rawCounts.whole.info[vids, ]

  }
  outdata$mrna.tissues <- tmp.mrna.tissues
  outdata$mrna.tissues.model2 <- tmp.mrna.tissues.model2
  outdata$mrna.tissues.adjRIN <- tmp.mrna.tissues.adjRIN
  outdata$mrna.tissues.org <- tmp.mrna.tissues.org
  outdata$mrna.tissues.info <- tmp.mrna.tissues.info
  return(outdata)
}

mRNA_batch_RIN_analysis <- function(rnaData, tissuenames, outpath, beta_cutoff = 0.008) {
  # use updated rnaData
  # plot  RIN and age
  tmpdata <- data.frame(Age = rnaData$mrna.whole.info$age,
    Agegroup = factor(rnaData$mrna.whole.info$type,
      levels = c("Juvenile", "Young_adult", "Middle_aged", "Elderly")),
    Counts = colSums(rnaData$rawCounts.whole)[rownames(rnaData$mrna.whole.info)],
    RIN = as.numeric(rnaData$mrna.whole.info$RIN))

  pdf(paste0(outpath, "/RIN_vs_age.pdf"), width = 5.5, height = 5)
  tmpcor <- cor.test(tmpdata$Age, tmpdata$RIN)
  pp <- ggplot(data = tmpdata, aes(x = Age, y = RIN, color = Age)) + theme_classic() +
    geom_jitter(width = 0.2) +
    ylab("RIN") + xlab("Age") +
    annotate("text", x = 13, y = 9.2,
      label = paste0("Rho = ", signif(tmpcor$estimate, 3), "; p = ", signif(tmpcor$p.value, 3)),
      size = 6, color = "black") +
    lghplot.addtheme(legend.position = "right", hjust = 1, size = 18) +
    ggtitle("RIN vs Age")
  # theme(legend.text=element_text(size=18))
  print(pp)
  dev.off()

  pdf(paste0(outpath, "/Counts_vs_age.pdf"), width = 6.5, height = 5)
  tmpcor <- cor.test(tmpdata$Age, tmpdata$Counts)
  pp <- ggplot(data = tmpdata, aes(x = Age, y = Counts, color = Age)) + theme_classic() +
    geom_jitter(width = 0.2) +
    ylab("Total counts") + xlab("Age") +
    annotate("text", x = 13, y = 2.7e7,
      label = paste0("Rho = ", signif(tmpcor$estimate, 3), "; p = ", signif(tmpcor$p.value, 3)),
      size = 6, color = "black") +
    lghplot.addtheme(legend.position = "right", hjust = 1, size = 18) +
    ggtitle("Total counts vs Age")
  # theme(legend.text=element_text(size=18))
  print(pp)
  dev.off()


  ## compare adjust and non adjust.

  tmplm.adjall <- get_tissue_DEgenes_lm(rnaData$mrna.tissues[tissuenames],
    rnaData$mrna.tissues.info[tissuenames], tissue.systems = NULL)

  tmplm.adjRIN <- get_tissue_DEgenes_lm(rnaData$mrna.tissues.adjRIN[tissuenames],
    rnaData$mrna.tissues.info[tissuenames], tissue.systems = NULL)
  tmplm.org <- get_tissue_DEgenes_lm(rnaData$mrna.tissues.org[tissuenames],
    rnaData$mrna.tissues.info[tissuenames], tissue.systems = NULL)
  tmplm.model2 <- get_tissue_DEgenes_lm(rnaData$mrna.tissues.model2[tissuenames],
    rnaData$mrna.tissues.info[tissuenames], tissue.systems = NULL)
  # compare adj or not adj

  pdf(paste0(outpath, "/adj_beta_corrleation.pdf"), width = 4.5, height = 4)
  vids <- intersect(tmplm.adjall$MetaLimma$ID, tmplm.org$MetaLimma$ID)
  tmprho <- cor.test(tmplm.adjall$MetaLimma[vids, ]$MetaBeta, tmplm.org$MetaLimma[vids, ]$MetaBeta)$estimate
  p1 <- ggplot(, aes(x = tmplm.adjall$MetaLimma[vids, ]$MetaBeta, y = tmplm.org$MetaLimma[vids, ]$MetaBeta)) +
    theme_classic() +
    geom_point(size = 0.2, alpha = 0.3) + geom_smooth(method = "lm", size = 0.5) + # xlim(c(-0.03,0.08)) +ylim(c(-0.04,0.08))+
    ylab("Age effect size (not adjust)") + xlab("Age effect size (adjust RIN & batch)") +
    annotate("text", x = -0.01, y = 0.05, label = paste0("Rho = ", signif(tmprho, 4)), size = 6, color = "black") +
    lghplot.addtheme(legend.position = "none", hjust = 1, size = 14) +
    ggtitle("Adjust RIN & batch vs no adjust")
  # theme(legend.text=element_text(size=18))
  print(p1)
  # dev.off()

  vids <- intersect(tmplm.adjRIN$MetaLimma$ID, tmplm.org$MetaLimma$ID)
  tmprho <- cor.test(tmplm.adjRIN$MetaLimma[vids, ]$MetaBeta, tmplm.org$MetaLimma[vids, ]$MetaBeta)$estimate
  p2 <- ggplot(, aes(x = tmplm.adjRIN$MetaLimma[vids, ]$MetaBeta, y = tmplm.org$MetaLimma[vids, ]$MetaBeta)) +
    theme_classic() +
    geom_point(size = 0.2, alpha = 0.3) + geom_smooth(method = "lm", size = 0.5) + # xlim(c(-0.03,0.08)) +ylim(c(-0.04,0.08))+
    ylab("Age effect size (not adjust)") + xlab("Age effect size (adjust RIN)") +
    annotate("text", x = -0.01, y = 0.05, label = paste0("Rho = ", signif(tmprho, 4)), size = 6, color = "black") +
    lghplot.addtheme(legend.position = "none", hjust = 1, size = 14) +
    ggtitle("Adjust RIN vs no adjust")
  # theme(legend.text=element_text(size=18))
  print(p2)


  vids <- intersect(tmplm.adjall$MetaLimma$ID, tmplm.adjRIN$MetaLimma$ID)
  tmprho <- cor.test(tmplm.adjall$MetaLimma[vids, ]$MetaBeta, tmplm.adjRIN$MetaLimma[vids, ]$MetaBeta)$estimate
  p3 <- ggplot(, aes(x = tmplm.adjall$MetaLimma[vids, ]$MetaBeta, y = tmplm.adjRIN$MetaLimma[vids, ]$MetaBeta)) +
    theme_classic() +
    geom_point(size = 0.2, alpha = 0.3) + geom_smooth(method = "lm", size = 0.5) + # xlim(c(-0.03,0.08)) +ylim(c(-0.04,0.08))+
    ylab("Age effect size (adjust RIN)") + xlab("Age effect size (adjust RIN & batch)") +
    annotate("text", x = -0.01, y = 0.05, label = paste0("Rho = ", signif(tmprho, 4)), size = 6, color = "black") +
    lghplot.addtheme(legend.position = "none", hjust = 1, size = 14) +
    ggtitle("Adjust RIN & batch vs Adjust RIN")
  # theme(legend.text=element_text(size=18))
  print(p3)
  dev.off()

  ## ggven up
  up.noadj <- tmplm.org$MetaLimma$ID[tmplm.org$MetaLimma$MetaFDR < 0.05 &
    tmplm.org$MetaLimma$MetaBeta > beta_cutoff &
    !tmplm.org$MetaLimma$manyNA]

  up.adjRIN <- tmplm.adjRIN$MetaLimma$ID[tmplm.adjRIN$MetaLimma$MetaFDR < 0.05 &
    tmplm.adjRIN$MetaLimma$MetaBeta > beta_cutoff &
    !tmplm.adjRIN$MetaLimma$manyNA]

  up.adjall <- tmplm.adjall$MetaLimma$ID[tmplm.adjall$MetaLimma$MetaFDR < 0.05 &
    tmplm.adjall$MetaLimma$MetaBeta > beta_cutoff &
    !tmplm.adjall$MetaLimma$manyNA]
  up.model2 <- tmplm.adjall$MetaLimma$ID[tmplm.model2$MetaLimma$MetaFDR < 0.05 &
    tmplm.model2$MetaLimma$MetaBeta > beta_cutoff &
    !tmplm.model2$MetaLimma$manyNA]

  pdf(paste0(outpath, "/ggven_identified_genes_up.pdf"), width = 7, height = 7)

  pp <- ggvenn::ggvenn(data = list(Up_noadj = up.noadj,
    Up_adjRIN = up.adjRIN,
    up_adjRIN_batch = up.adjall),
  fill_color = c("cornflowerblue", "green", "#E64B35FF"),
  show_percentage = F, stroke_size = 0.5, stroke_alpha = 0.6, text_size = 9
  )
  print(pp)
  dev.off()

  pdf(paste0(outpath, "/ggven_identified_genes_up_model1_vs_model2.pdf"), width = 5, height = 3.5)

  pp <- ggvenn::ggvenn(data = list(Up_model1 = up.adjall,
    Up_model2 = up.model2),
  fill_color = c("#E64B35FF", "cornflowerblue"),
  show_percentage = F, stroke_size = 0.5, stroke_alpha = 0.6, text_size = 9
  )
  print(pp)
  dev.off()



  ## ggven down
  down.noadj <- tmplm.org$MetaLimma$ID[tmplm.org$MetaLimma$MetaFDR < 0.05 &
    tmplm.org$MetaLimma$MetaBeta < -beta_cutoff &
    !tmplm.org$MetaLimma$manyNA]

  down.adjRIN <- tmplm.adjRIN$MetaLimma$ID[tmplm.adjRIN$MetaLimma$MetaFDR < 0.05 &
    tmplm.adjRIN$MetaLimma$MetaBeta < -beta_cutoff &
    !tmplm.adjRIN$MetaLimma$manyNA]

  down.adjall <- tmplm.adjall$MetaLimma$ID[tmplm.adjall$MetaLimma$MetaFDR < 0.05 &
    tmplm.adjall$MetaLimma$MetaBeta < -beta_cutoff &
    !tmplm.adjall$MetaLimma$manyNA]
  down.model2 <- tmplm.adjall$MetaLimma$ID[tmplm.model2$MetaLimma$MetaFDR < 0.05 &
    tmplm.model2$MetaLimma$MetaBeta < -beta_cutoff &
    !tmplm.model2$MetaLimma$manyNA]

  pdf(paste0(outpath, "/ggven_identified_genes_down.pdf"), width = 7, height = 7)

  pp <- ggvenn::ggvenn(data = list(Down_noadj = down.noadj,
    Down_adjRIN = down.adjRIN,
    Down_adjRIN_batch = down.adjall),
  fill_color = c("cornflowerblue", "green", "#E64B35FF"),
  show_percentage = F, stroke_size = 0.5, stroke_alpha = 0.6, text_size = 9
  )
  print(pp)
  dev.off()

  pdf(paste0(outpath, "/ggven_identified_genes_down_model1_vs_model2.pdf"), width = 5, height = 3)

  pp <- ggvenn::ggvenn(data = list(Down_model1 = down.adjall,
    Down_model2 = down.model2),
  fill_color = c("#E64B35FF", "cornflowerblue"),
  show_percentage = F, stroke_size = 0.5, stroke_alpha = 0.6, text_size = 9
  )
  print(pp)
  dev.off()


  pdf(paste0(outpath, "/ggven_identified_genes_model1_vs_model2.pdf"), width = 7, height = 7)

  pp <- ggvenn::ggvenn(data = list(Up_model1 = up.adjall,
    Up_model2 = up.model2,
    Down_model1 = down.adjall,
    Down_model2 = down.model2),
  fill_color =  c("darkorchid1", "yellow", "green", "cornflowerblue"),
  show_percentage = F, stroke_size = 0.5, stroke_alpha = 0.6, text_size = 9
  )
  print(pp)
  dev.off()
  return(1)
}

updateRNA_withLimma <- function(rnaData, headers, tissue.systems, macacageneInfo, outpath) {
  # filter with Rin > 6 and mapped counts > 10M
  idx <- colSums(rnaData$rawCounts.whole) >= 10000000 & rnaData$rawCounts.whole.info$RIN >= 6
  rnaData$rawCounts.whole <- rnaData$rawCounts.whole[, idx]
  rnaData$rawCounts.whole.info <- rnaData$rawCounts.whole.info[idx, ]

  outdata <- rnaData
  tcounts <- delete_dup_genes(rnaData$rawCounts.whole, headers)
  tcounts <- tcounts[rownames(tcounts) != "", ]

  # consider coding genes
  tmptype <- macacageneInfo[rownames(tcounts), ]$Type
  idx <- tmptype == "protein_coding" #| substr(tmptype,nchar(tmptype)-4, nchar(tmptype)) == "_gene"
  idx[is.na(idx)] <- FALSE
  tcounts <- tcounts[idx, ]
  tcpm <- get_CPM_counts(tcounts)
  tcpm.nonormal <- get_CPM_nonormal(tcounts)

  batch <- rnaData$rawCounts.whole.info$batch
  RINs <- as.numeric(rnaData$rawCounts.whole.info$RIN)
  RINs[is.na(RINs)] <- mean(RINs, na.rm = T)

  #
  metadata <- data.frame(batch = batch,
    RIN = RINs,
    age = rnaData$rawCounts.whole.info$age,
    tissue = rnaData$rawCounts.whole.info$tissues)
  metadata$batch <- factor(metadata$batch)
  metadata$tissue <- factor(metadata$tissue)


  # Define the batch variable
  batch <- metadata$batch

  tcpm.batchadj <- limma::removeBatchEffect(tcpm, batch = batch, covariates = metadata$RIN)



  # for tsne plot Figure 1
  mrna.whole.std <- standardise_matrix(tcpm)
  mrna.whole.std.adj <- standardise_matrix(tcpm.batchadj)
  # mrna.whole.std.adj1 = standardise_matrix(tcpm.batchadj1)

  ## PCA before and after
  # because the expression data  is locate in the same magnitude,
  pca_before <- prcomp(t(tcpm), cor = F)
  percPCA_before <- 100 * summary(pca_before)$importance

  pca_after <- prcomp(t(tcpm.batchadj), cor = F)
  percPCA_after <- 100 * summary(pca_after)$importance

  ## boxplot
  png(paste0(outpath, "/boxplot_mRNA_before_adj.png"), width = 2000, height = 500)
  boxplot(tcpm, ylab = "Normalized log2 CPM", cex = 0.2,
    xlab = paste0("Samples of all mRNAs", " (n= ", ncol(tcpm), ")")
  )
  dev.off()

  ## boxplot
  png(paste0(outpath, "/boxplot_mRNA_after_adj.png"), width = 2000, height = 500)
  boxplot(tcpm.batchadj, ylab = "Normalized log2 CPM", cex = 0.2,
    xlab = paste0("Samples of all mRNAs after batch adj", " (n= ", ncol(tcpm.batchadj), ")")
  )
  dev.off()


  # plot PC1 PC2 and batch
  pdf(paste0(outpath, "/CPM_PC1_vs_PC2_before_batch_adj.pdf"), width = 9, height = 7)
  pp <- ggplot(, aes(pca_before$x[, 1], pca_before$x[, 2], color = batch)) + geom_point() +
    theme_bw() + lghplot.addtheme(legend.position = "right") +
    xlab(paste("PC1(", as.character(round(percPCA_before[2, 1], 1)), "%)", sep = "")) +
    ylab(paste("PC2(", as.character(round(percPCA_before[2, 2], 1)), "%)", sep = "")) +
    labs(title = paste0("PC1 vs PC2 before batch adjusted in mRNA"))
  print(pp)
  dev.off()

  # plot PC1 PC2  after batch
  pdf(paste0(outpath, "/CPM_PC1_vs_PC2_after_batch_adj.pdf"), width = 9, height = 7)
  pp <- ggplot(, aes(pca_after$x[, 1], pca_after$x[, 2], color = batch)) + geom_point() +
    theme_bw() + lghplot.addtheme(legend.position = "right") +
    xlab(paste("PC1(", as.character(round(percPCA_after[2, 1], 1)), "%)", sep = "")) +
    ylab(paste("PC2(", as.character(round(percPCA_after[2, 2], 1)), "%)", sep = "")) +
    labs(title = paste0("PC1 vs PC2 after batch adjusted in mRNA"))
  print(pp)
  dev.off()


  ## plot RIN vs PC1 or PC2
  pdf(paste0(outpath, "/PC1_vs_RIN_before_adj.pdf"), width = 9, height = 7)
  tmpxx <- cor.test(pca_before$x[, 1], RINs)
  tmptext <- paste0("r=", signif(tmpxx$estimate, 2), "; p=", signif(tmpxx$p.value, 2))
  pp <- ggplot(, aes(pca_before$x[, 1], RINs, color = batch)) +
    geom_point() + scale_color_npg() +
    theme_bw() + lghplot.addtheme(legend.position = "right") + # stat_ellipse(lwd=1,level = 0.99) +
    annotate(geom = "text", x = summary(pca_before$x[, 1])[5],
      y = summary(RINs)[5], label = tmptext,
      color = "black", size = 6) +
    xlab(paste("PC1(", as.character(round(percPCA_before[2, 1], 1)), "%)", sep = "")) +
    ylab("RIN value")
  print(pp)
  dev.off()

  pdf(paste0(outpath, "/PC2_vs_RIN_before_adj.pdf"), width = 9, height = 7)
  tmpxx <- cor.test(pca_before$x[, 2], RINs)
  tmptext <- paste0("r=", signif(tmpxx$estimate, 2), "; p=", signif(tmpxx$p.value, 2))
  pp <- ggplot(, aes(pca_before$x[, 2], RINs, color = batch)) +
    geom_point() + scale_color_npg() +
    theme_bw() + lghplot.addtheme(legend.position = "right") +
    annotate(geom = "text", x = summary(pca_before$x[, 2])[5],
      y = summary(RINs)[5], label = tmptext,
      color = "black", size = 6) +
    xlab(paste("PC2(", as.character(round(percPCA_before[2, 2], 1)), "%)", sep = "")) +
    ylab("RIN value")
  print(pp)
  dev.off()

  ## plot RIN vs PC1 or PC2 after batch adj
  pdf(paste0(outpath, "/PC1_vs_RIN_after_adj.pdf"), width = 9, height = 7)
  tmpxx <- cor.test(pca_after$x[, 1], RINs)
  tmptext <- paste0("r=", signif(tmpxx$estimate, 2), "; p=", signif(tmpxx$p.value, 2))
  pp <- ggplot(, aes(pca_after$x[, 1], RINs, color = batch)) +
    geom_point() + scale_color_npg() +
    theme_bw() + lghplot.addtheme(legend.position = "right") +
    annotate(geom = "text", x = summary(pca_after$x[, 1])[5],
      y = summary(RINs)[5], label = tmptext,
      color = "black", size = 6) +
    xlab(paste("PC1(", as.character(round(percPCA_after[2, 1], 1)), "%)", sep = "")) +
    ylab("RIN value")
  print(pp)
  dev.off()

  pdf(paste0(outpath, "/PC2_vs_RIN_after_adj.pdf"), width = 9, height = 7)
  tmpxx <- cor.test(pca_after$x[, 2], RINs)
  tmptext <- paste0("r=", signif(tmpxx$estimate, 2), "; p=", signif(tmpxx$p.value, 2))
  pp <- ggplot(, aes(pca_after$x[, 2], RINs, color = batch)) +
    geom_point() + scale_color_npg() +
    theme_bw() + lghplot.addtheme(legend.position = "right") +
    annotate(geom = "text", x = summary(pca_after$x[, 2])[5],
      y = summary(RINs)[5], label = tmptext,
      color = "black", size = 6) +
    xlab(paste("PC2(", as.character(round(percPCA_after[2, 2], 1)), "%)", sep = "")) +
    ylab("RIN value")
  print(pp)
  dev.off()

  #
  ### for tsne for batch before adjust
  p1 <- tsne(mrna.whole.std, labels = tissue.systems[rnaData$rawCounts.whole.info$tissues],
    legendtextsize = 10, dotsize = 2, seed = 2025520)
  p1 <- p1 + theme_classic() + lghplot.addtheme(legend.position = "none") +
    scale_color_npg() + xlab("tsne 1") + ylab("tsne 2") +
    theme(axis.line = element_line(size = 1.0)) + ggtitle("Transcriptome")
  pdf(file = paste0(outpath, "/tSNE_mrna_using_standardised_before_adjust.pdf"), height = 4, width = 4)
  print(p1)
  dev.off()


  ### for tsne for batch after ajust
  p2 <- tsne(mrna.whole.std.adj, labels = tissue.systems[rnaData$rawCounts.whole.info$tissues],
    legendtextsize = 10, dotsize = 2, seed = 2025520)
  p2 <- p2 + theme_classic() + lghplot.addtheme(legend.position = "none") +
    scale_color_npg() + xlab("tsne 1") + ylab("tsne 2") +
    theme(axis.line = element_line(size = 1.0)) + ggtitle("Transcriptome")
  pdf(file = paste0(outpath, "/tSNE_mrna_using_standardised_after_adjust.pdf"), height = 4, width = 4)
  print(p2)
  dev.off()

  outdata$mrna.whole <- tcpm.batchadj
  outdata$mrna.whole.info <- outdata$rawCounts.whole.info

  # tissues
  tmp.mrna.tissues <- list()
  tmp.mrna.tissues.info <- list()
  tmp.tissue.names <- unique(outdata$rawCounts.whole.info$tissues)
  for (i in 1:length(tmp.tissue.names)) {
    thisname <- tmp.tissue.names[i]
    vids <- outdata$rawCounts.whole.info$tissues == thisname
    # get high confidance tissue data
    thiscouts <- tcounts[, vids]
    # thisinfo = outdata$rawCounts.whole.info[vids,]
    tmpcpm.nonormal <- tcpm.nonormal[, vids]
    xids <- rowMeans(tmpcpm.nonormal) > 1 & rowSums(thiscouts < 5) < ncol(thiscouts) * 0.2
    tmp.mrna.tissues[[thisname]] <- tcpm.batchadj[xids, vids]
    tmp.mrna.tissues.info[[thisname]] <- outdata$rawCounts.whole.info[vids, ]
    tmpthecounts <- outdata$rawCounts.whole.info[vids, ]
  }
  outdata$mrna.tissues <- tmp.mrna.tissues
  outdata$mrna.tissues.info <- tmp.mrna.tissues.info
  return(outdata)
}

## meta
short_met_names <- function(Metamet.typeI, outype = "data") {
  tmpmet <- Metamet.typeI
  sortname <- rownames(tmpmet)
  sortname <- gsub("N-Acetyl-α-D-glucosamine 1-phosphate", "GlcNAc-1-phosphate", sortname)
  sortname <- gsub("Flavin mononucleotide \\(FMN\\)", "Flavin mononucleotide", sortname)
  sortname <- gsub("8Z,11Z,14Z-Eicosatrienoic acid", "Eicosatrienoic acid", sortname)
  sortname <- gsub("3-phenethyl-2-thioxoimidazolidin-4-one", "Pubchem:2726668", sortname)
  sortname <- gsub("2-\\(1,3-dimethyl-1H-pyrazol-5-yl\\)-1H-isoindole-1,3\\(2H\\)-dione", "CHEMBL1705433", sortname)
  sortname <- gsub("1-\\[4-\\(1-adamantyl\\)phenoxy\\]-3-piperidinopropan-2-ol\\s+hydrochloride", "CAS175136-32-0", sortname)
  sortname <- gsub("2-\\[4-\\(tert-butyl\\)-2,3-dihydro-1,3-thiazol-2-yliden]-3-oxobutanenitrile", "CB81536810", sortname)
  sortname <- gsub("3,4-Dihydroxyphenylpropionic acid", "Dihydrocaffeic acid", sortname)
  sortname <- gsub("ethyl 4-methyl-2-\\(1H-pyrrol-1-yl\\)-1,3-thiazole-5-carboxylate", "CHEBI:190645", sortname)
  sortname[sortname == "δ-Valerolactam"] <- "Deta-Valerolactame"
  sortname[sortname == "4,4'-dimethoxy[1,1'-biphenyl]-2-carbonitrile"] <- "CHEBI:183677"
  sortname[sortname == "4-(3-Hydroxybutyl)phenyl β-D-glucopyranoside"] <- "CHEBI:167790"
  sortname <- gsub("Acetyl-\\β-methylcholine", "Acetyl-beta-methylcholine", sortname)
  sortname[sortname == "3-(methylsulfonyl)-2H-chromen-2-one"] <- "CHEBI:192673"
  sortname[sortname == "2-[4-(tert-butyl)-2,3-dihydro-1,3-thiazol-2-yliden]-3-oxobutanenitrile"] <- "CAS:307975-76-4"
  sortname[sortname == "3-phenethyl-2-thioxoimidazolidin-4-one"] <- "CAS:37021-14-0"
  sortname[sortname == "2-{2-[4-(tert-butyl)phenyl]-4-hydroxy-1,3-thiazol-5-yl}acetic acid"] <- "CHEBI:183477"
  sortname[sortname == "5-fluoro AB-PINACA N-(4-hydroxypentyl) metabolite"] <- "CAS:2460433-23-0"
  sortname[sortname == "3-(2,3-dihydro-1H-indol-1-yl)-2-[(2-furylmethyl)sulfonyl]acrylonitrile"] <- "Alpha-methyldeoxybenzoin"
  sortname[sortname == "N-[(4-hydroxy-3-methoxyphenyl)methyl]-8-methylnonanamide"] <- "CHEMBL419862"
  sortname[sortname == "2-{[2-(4-methylpiperazino)phenyl]methylene}hydrazine-1-carbothioamide"] <- "CAS:2929-81-9"
  sortname[sortname == "6-methyl-4-(morpholinomethyl)-2H-chromen-2-one"] <- "CHEBI:190578"
  sortname[sortname == "3-(2-methylpropyl)-octahydropyrrolo[1,2-a]pyrazine-1,4-dione"] <- "CHEBI:228605"
  sortname[sortname == "3-(4-Hydroxyphenyl)propionic acid"] <- "Phloretic acid"
  sortname[sortname == "16α-Hydroxydehydroepiandrosterone"] <- "16alpha-OH-DHEA"
  sortname[sortname == "11-α-Hydroxy-17-methyltestosterone"] <- "CHEBI:34138"
  sortname[sortname == "7α-Hydroxytestosterone"] <- "7a-Hydroxytestosterone"
  sortname[sortname == "2-Methylbutyl beta-D-glucopyranoside"] <- "beta-Methylglucoside"
  sortname[sortname == "ethyl 4-methyl-2-(1H-pyrrol-1-yl)-1,3-thiazole-5-carboxylate"] <- "CHEBI:190645"
  sortname[sortname == "3'-Adenosine monophosphate (3'-AMP)"] <- "3'-AMP"
  sortname[sortname == "N-Acetyl-1-aspartylglutamic acid"] <- "NAAG"
  sortname[sortname == "α-Hydroxyhippuric acid"] <- "a-Hydroxyhippuric acid"
  sortname[sortname == "Phosphatidylinositol-1,2-dipalmitoyl"] <- "PI(16:0/16:0)"
  sortname[sortname == "Cytidine-5'-monophosphate"] <- "CMP"
  sortname[sortname == "alpha-D-Glucopyranosyl 2-O-(2-methylbutanoyl)-alpha-D-glucopyranoside"] <- "CHEBI:183682"
  sortname[sortname == "2-(2-chlorophenyl)-1-cyclohexyl-6-oxopiperidine-3-carboxylic acid"] <- "CAS:1212241-27-4"
  sortname[sortname == "Adenosine 5'-monophosphate"] <- "AMP"
  sortname[sortname == "XLR11 N-(2-fluoropentyl) isomer"] <- "CHEBI:190582"
  sortname[sortname == "2-(tert-butyl)-1,3-thiazolane-4-carboxylic acid"] <- "CAS:1012881-39-8"
  sortname[sortname == "XLR11 N-(2-fluoropentyl) isomer"] <- "CHEBI:190582"
  sortname[sortname == "2-(acetylamino)-3-(1H-indol-3-yl)propanoic acid"] <- "CAS:2280-01-5"
  sortname[sortname == "α-Ethylaminopentiophenone"] <- "a-Ethylaminopentiophenone"
  rownames(tmpmet) <- sortname
  if (tolower(outype) == "data") {
    return(tmpmet)
  } else {
    return(sortname)
  }
}


lmGenes <- function(expr, factor, cutn = 5) {
  ## for age
  require(limma)
  design <- model.matrix(~factor)
  fit <- lmFit(expr, design)
  fit <- eBayes(fit)
  results <- topTable(fit, adjust.method = "BH", number = Inf)
  tmpu <- function(x) {
    length(unique(x))
  }
  idna <- apply(expr[rownames(results), ], 1, tmpu) < cutn
  results[idna, ] <- NA
  results  <- cbind(data.frame(ID = rownames(results)), results)
  colnames(results)[colnames(results) == "logFC"] <- "beta"
  colnames(results)[colnames(results) == "P.Value"] <- "Pvalue"
  colnames(results)[colnames(results) == "adj.P.Val"] <- "FDR"
  return(results)
}


get_tissue_DEgenes_lm <- function(tissues.expr, tissues.clin, tissue.systems = tissue.systems) {
  # 1.  tissue aging related genes using lm expr~age
  print("Step 1: get age related mols by lm expr~age for each tissue")
  alltissues <- names(tissues.expr)

  DEprobeta <- list()
  DEproPvalue <- list()
  DEproAging <- list()
  for (i in 1:length(alltissues)) {
    thistissue <- alltissues[i]
    thispro <- tissues.expr[[thistissue]]
    thispro.info <- tissues.clin[[thistissue]]
    thisDEpro <- lmGenes(thispro, thispro.info$age, cutn = ncol(thispro) / 2)
    DEproAging[[i]] <- thisDEpro
    DEprobeta[[i]] <- thisDEpro$beta
    names(DEprobeta[[i]]) <- rownames(thisDEpro)
    DEproPvalue[[i]] <- thisDEpro$Pvalue
    names(DEproPvalue[[i]])  <- rownames(thisDEpro)
  }
  names(DEproAging) <- alltissues
  names(DEprobeta) <- alltissues
  names(DEproPvalue) <- alltissues

  # num up and down proteins
  DEprobeta_matrix <- list_to_matrix(DEprobeta, alltissues)
  DEproPvalue_matrix <- list_to_matrix(DEproPvalue, alltissues)

  Aging_pro_sigup_matrix <- (DEprobeta_matrix > 0 & DEproPvalue_matrix < 0.05) + 0
  Aging_pro_sigdown_matrix <- -((DEprobeta_matrix < -0 & DEproPvalue_matrix < 0.05) + 0)
  Aging_pro_sigall_matrix <- Aging_pro_sigup_matrix + Aging_pro_sigdown_matrix

  if (!is.null(tissue.systems)) {
    Aging_pro_updown <- data.frame(stringsAsFactors = F,
      num.up = colSums(Aging_pro_sigup_matrix, na.rm = T) /
        colSums(!is.na(Aging_pro_sigup_matrix)),
      num.down = colSums(Aging_pro_sigdown_matrix, na.rm = T) /
        colSums(!is.na(Aging_pro_sigdown_matrix)),
      num.all = colSums(abs(Aging_pro_sigall_matrix), na.rm = T) /
        colSums(!is.na(Aging_pro_sigall_matrix)),
      tissues = colnames(Aging_pro_sigup_matrix),
      tissue_systems = tissue.systems)
    rownames(Aging_pro_updown) <- Aging_pro_updown$tissues

  } else {

    Aging_pro_updown <- data.frame(stringsAsFactors = F,
      num.up = colSums(Aging_pro_sigup_matrix, na.rm = T) /
        colSums(!is.na(Aging_pro_sigup_matrix)),
      num.down = colSums(Aging_pro_sigdown_matrix, na.rm = T) /
        colSums(!is.na(Aging_pro_sigdown_matrix)),
      num.all = colSums(abs(Aging_pro_sigall_matrix), na.rm = T) /
        colSums(!is.na(Aging_pro_sigall_matrix)),
      tissues = colnames(Aging_pro_sigup_matrix))
    rownames(Aging_pro_updown) <- Aging_pro_updown$tissues

  }


  outlist <- list()
  outlist$Aging <- DEproAging
  outlist$beta <- DEprobeta
  outlist$Pvalue <- DEproPvalue
  outlist$Aging_updown <- Aging_pro_updown



  # 2. set data for common age related genes across tissues
  print("Step 2: construct data for meta analysis")
  exprData <- list()
  allgenes <- c()
  thenames <- names(tissues.clin)

  for (i in 1:length(thenames)) {
    studyname <- thenames[i]
    # expr.log2 = #log2(repExpr[[studyname]]+1e-6)
    expr.log2 <- tissues.expr[[studyname]]
    allgenes <- unique(c(allgenes, rownames(expr.log2)))
    exprData[[i]] <- expr.log2
  }
  names(exprData) <- thenames
  # filter allgenes NA
  nas  <- rep(0, length(allgenes))
  k <- 0
  for (i in 1:length(thenames)) {
    tmp <- as.data.frame(exprData[[i]])
    tmp <- as.matrix(tmp[allgenes, ])
    k <- k + ncol(tmp)
    nas <- nas + rowSums(is.na(tmp))
  }
  nas <- nas / k
  names(nas) <- allgenes

  data.AveExpr <- matrix(0, length(allgenes), length(tissues.clin))
  data.beta <- matrix(0, length(allgenes), length(tissues.clin))
  data.Pvalue <- matrix(0, length(allgenes), length(tissues.clin))
  data.FDR <- matrix(0, length(allgenes), length(tissues.clin))
  for (i in 1:length(thenames)) {
    studyname <- thenames[i]
    tmp <- as.data.frame(exprData[[i]])
    tmp <- as.matrix(tmp[allgenes, ])
    tmp1 <- fillgaps_rowMinNA_1(tmp)
    rownames(tmp1) <- allgenes

    tmplm <- outlist$Aging[[i]]
    tmplm <- tmplm[allgenes, ]

    exprData[[i]] <- tmp1

    tmpclin <- tissues.clin[[i]]
    thisclin <- data.frame(age = tmpclin$age,
      type = tmpclin$type,
      tissue = rep(studyname, nrow(tmpclin)))

    if (i == 1) {
      exprData.matrix <- tmp1
      combine_clin <- thisclin
    } else {
      exprData.matrix <- cbind(exprData.matrix, tmp1)
      combine_clin <- rbind(combine_clin, thisclin)
    }

    data.beta[, i] <- tmplm$beta
    data.Pvalue[, i] <- tmplm$Pvalue
    data.FDR[, i] <- tmplm$FDR
    data.AveExpr[, i]  <-  tmplm$AveExpr
  }
  colnames(data.AveExpr) <- paste0("AveExpr_", thenames)
  colnames(data.beta) <- paste0("beta_", thenames)
  colnames(data.Pvalue) <- paste0("Pvalue_", thenames)
  colnames(data.FDR) <- paste0("FDR_", thenames)

  # 3. perform lm to get age-related genes across tissues use limma linear model
  #   combined_expr ~ age + tissue
  # Example: Differential Expression Analysis using limma
  print("Step 3: meta analysis using limma: combined_expr ~ age + tissue")

  design <- model.matrix(~ age + tissue, data = combine_clin)
  fit <- lmFit(exprData.matrix, design)
  fit <- eBayes(fit)
  results_age <- topTable(fit, coef = "age", adjust.method = "BH", number = Inf)
  results_age <- results_age[allgenes, ]

  MetaLimma <- data.frame(ID = allgenes, stringsAsFactors = F,
    MetaAveExpr = results_age$AveExpr,
    MetaBeta = results_age$logFC,
    MetaPvalue = results_age$P.Value,
    MetaFDR = results_age$adj.P.Val,
    row.names = allgenes)
  MetaLimma$MetaFDR[is.na(MetaLimma$MetaFDR)] <- 1
  MetaLimma$MetaPvalue[is.na(MetaLimma$MetaPvalue)] <- 1


  MetaLimma$NA_perc <- nas
  MetaLimma$manyNA <- rowSums(is.na(data.beta)) >= ncol(data.beta) / 3
  MetaLimma <- cbind(MetaLimma, data.beta)
  MetaLimma <- cbind(MetaLimma, data.Pvalue)
  MetaLimma <- cbind(MetaLimma, data.FDR)

  outlist$MetaLimma <- MetaLimma
  outlist$LimmaFit <- fit

  print("Finished!")

  return(outlist)
}


obtain_updown_mols <- function(DEpro.tissues.lm, betacutoff = 0.008, tissue.systems = tissue.systems) {

  DEprobeta  <- DEpro.tissues.lm$beta
  DEproPvalue  <- DEpro.tissues.lm$Pvalue
  alltissues <- names(DEpro.tissues.lm$beta)
  # num up and down proteins
  DEprobeta_matrix <- list_to_matrix(DEprobeta, alltissues)
  DEproPvalue_matrix <- list_to_matrix(DEproPvalue, alltissues)

  Aging_pro_sigup_matrix <- (DEprobeta_matrix > betacutoff & DEproPvalue_matrix < 0.05) + 0
  Aging_pro_sigdown_matrix <- -((DEprobeta_matrix < -betacutoff & DEproPvalue_matrix < 0.05) + 0)
  Aging_pro_sigall_matrix <- Aging_pro_sigup_matrix + Aging_pro_sigdown_matrix

  if (!is.null(tissue.systems)) {
    Aging_pro_updown <- data.frame(stringsAsFactors = F,
      num.up = colSums(Aging_pro_sigup_matrix, na.rm = T) /
        colSums(!is.na(Aging_pro_sigup_matrix)),
      num.down = colSums(Aging_pro_sigdown_matrix, na.rm = T) /
        colSums(!is.na(Aging_pro_sigdown_matrix)),
      num.all = colSums(abs(Aging_pro_sigall_matrix), na.rm = T) /
        colSums(!is.na(Aging_pro_sigall_matrix)),
      tissues = colnames(Aging_pro_sigup_matrix),
      tissue_systems = tissue.systems)
    rownames(Aging_pro_updown) <- Aging_pro_updown$tissues

  } else {

    Aging_pro_updown <- data.frame(stringsAsFactors = F,
      num.up = colSums(Aging_pro_sigup_matrix, na.rm = T) /
        colSums(!is.na(Aging_pro_sigup_matrix)),
      num.down = colSums(Aging_pro_sigdown_matrix, na.rm = T) /
        colSums(!is.na(Aging_pro_sigdown_matrix)),
      num.all = colSums(abs(Aging_pro_sigall_matrix), na.rm = T) /
        colSums(!is.na(Aging_pro_sigall_matrix)),
      tissues = colnames(Aging_pro_sigup_matrix))
    rownames(Aging_pro_updown) <- Aging_pro_updown$tissues

  }
  return(Aging_pro_updown)

}



plot_save_DEgenes <- function(DEpro.tissues.lm,
                              pdffile,
                              pngfile,
                              excelfile, cutoff = 0.008) {
  proteinVoconoPlot <- list()
  names.DEproAging <- names(DEpro.tissues.lm$Aging)
  for (i in 1:length(DEpro.tissues.lm$Aging)) {
    res2 <- DEpro.tissues.lm$Aging[[i]]
    res2 <- short_met_names(res2)
    res2$ID <- rownames(res2)

    proteinVoconoPlot[[i]] <- plot_DEflux(res2, title = names.DEproAging[i],
      x = "beta", y = "Pvalue",
      xlab = expression("Age effect size (" * beta[t] * ")"), 
      ylab = bquote(~ -log[10] ~ italic(P)),
      FCcutoff = cutoff, num.showlab = 3)
  }

  pdf(file = pdffile, width = 3 * 5 + 1, height = 3 * 6)
  grid.arrange(arrangeGrob(grobs = proteinVoconoPlot, ncol = 5))
  dev.off()

  sheetNames <- names(DEpro.tissues.lm$Aging)
  xx <- list()
  for (i in 1:length(sheetNames)) {
    tmp <- DEpro.tissues.lm$Aging[[i]]
    tmp <- tmp[order(tmp$Pvalue), ]

    tmp[c(2, 3, 5)] <- signif(tmp[c(2, 3, 5)], 4)
    tmp <- tmp[, c(1, 2, 3, 5)]

    xx[[i]] <- tmp
  }
  names(xx) <- sheetNames
  openxlsx::write.xlsx(xx, file = excelfile, keepNA = T)
  return(1)
}

ggvenn_two_meta_methods <- function(DEpro.tissues.lm, outfile) {
  tmpaa <- DEpro.tissues.lm$MetaLimma

  pdf(outfile, width = 8)
  pp <- ggvenn::ggvenn(data = list(Up_Limma = tmpaa$ID[tmpaa$MetaFDR < 0.05 & tmpaa$MetaBeta > 0 & tmpaa$manyNA == FALSE],
    Up_MetaDE = tmpaa$ID[tmpaa$MetaFDR1 < 0.05 & tmpaa$MetaBeta1 > 0 & tmpaa$manyNA == FALSE],
    Down_Limma = tmpaa$ID[tmpaa$MetaFDR < 0.05 & tmpaa$MetaBeta < 0 & tmpaa$manyNA == FALSE],
    Down_MetaDE = tmpaa$ID[tmpaa$MetaFDR1 < 0.05 & tmpaa$MetaBeta1 < 0 & tmpaa$manyNA == FALSE]),
  fill_color = c("darkorchid1", "yellow", "green", "cornflowerblue"),
  show_percentage = F, stroke_size = 0.5, stroke_alpha = 0.6, text_size = 9)
  print(pp)
  dev.off()

  out <- list()
  out$commonup005 <- intersect(tmpaa$ID[tmpaa$MetaFDR < 0.05 & tmpaa$MetaBeta > 0 & tmpaa$manyNA == FALSE],
    tmpaa$ID[tmpaa$MetaFDR1 < 0.05 & tmpaa$MetaBeta1 > 0 & tmpaa$manyNA == FALSE])
  out$commondown005 <- intersect(tmpaa$ID[tmpaa$MetaFDR < 0.05 & tmpaa$MetaBeta < 0 & tmpaa$manyNA == FALSE],
    tmpaa$ID[tmpaa$MetaFDR1 < 0.05 & tmpaa$MetaBeta1 < 0 & tmpaa$manyNA == FALSE])
  out$commonup001 <- intersect(tmpaa$ID[tmpaa$MetaFDR < 0.01 & tmpaa$MetaBeta > 0 & tmpaa$manyNA == FALSE],
    tmpaa$ID[tmpaa$MetaFDR1 < 0.01 & tmpaa$MetaBeta1 > 0 & tmpaa$manyNA == FALSE])
  out$commondown001 <- intersect(tmpaa$ID[tmpaa$MetaFDR < 0.01 & tmpaa$MetaBeta < 0 & tmpaa$manyNA == FALSE],
    tmpaa$ID[tmpaa$MetaFDR1 < 0.01 & tmpaa$MetaBeta1 < 0 & tmpaa$manyNA == FALSE])

  return(out)
}


metaGene <- function(repExpr, repClin,
                     meta.method = "REM", tail = "abs", parametric = TRUE) {
  set.seed(123) #for comparison with previous
  
  #subfunc
  sub_fillgaps_MinNA_formeta <- function(X,nas,nascutoff = 0.3){
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
  
  sub_fillgaps_rowMinNA <- function(X){
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
  
  # set data
  flux3.log2 <- list()

  allrxns <- c()
  thenames <- names(repClin)

  for (i in 1:length(thenames)) {
    studyname <- thenames[i]
    # expr.log2 = #log2(repExpr[[studyname]]+1e-6)
    expr.log2 <- repExpr[[studyname]]
    allrxns <- unique(c(allrxns, rownames(expr.log2)))
    flux3.log2[[i]] <- expr.log2
  }
  # filter allrxns NA
  nas  <- rep(0, length(allrxns))
  k <- 0
  for (i in 1:length(thenames)) {
    tmp <- as.data.frame(flux3.log2[[i]])
    tmp <- as.matrix(tmp[allrxns, ])
    k <- k + ncol(tmp)
    nas <- nas + rowSums(is.na(tmp))
  }
  nas <- nas / k
  names(nas) <- allrxns

  # allrxns = allrxns[nas/k < 0.5]

  flux3.AveExpr <- matrix(0, length(allrxns), length(repClin))
  flux3.log2FC <- matrix(0, length(allrxns), length(repClin))
  flux3.Pvalue <- matrix(0, length(allrxns), length(repClin))
  flux3.FDR <- matrix(0, length(allrxns), length(repClin))
  for (i in 1:length(thenames)) {
    tmp <- as.data.frame(flux3.log2[[i]])
    tmp <- as.matrix(tmp[allrxns, ])
    tmp1 <- sub_fillgaps_rowMinNA(tmp)
    tmp <- sub_fillgaps_MinNA_formeta(tmp, nas)
    # tmp[is.na(tmp)] = 0
    rownames(tmp) <- allrxns
    rownames(tmp1) <- allrxns

    tmpDEflux <- DEGenes.simplified(tmp1, catagory = repClin[[i]]$type == "Elderly",
      subset = repClin[[i]]$type == "Elderly" | repClin[[i]]$type == "Young_adult")
    flux3.log2[[i]] <- tmp
    flux3.log2FC[, i] <- tmpDEflux$log2FC
    flux3.Pvalue[, i] <- tmpDEflux$Pvalue
    flux3.FDR[, i] <- tmpDEflux$FDR
    flux3.AveExpr[, i]  <-  tmpDEflux$AveExpr
  }
  colnames(flux3.AveExpr) <- paste0("AveExpr_", thenames)
  colnames(flux3.log2FC) <- paste0("log2FC_", thenames)
  colnames(flux3.Pvalue) <- paste0("Pvalue_", thenames)
  colnames(flux3.FDR) <- paste0("FDR_", thenames)

  # run metaDE
  data <- flux3.log2
  label <- list()
  K <- length(data)
  for (i in 1:length(repClin)) {
    label[[i]] <- repClin[[i]]$type
  }
  clin.data <- lapply(label, function(x) {
    data.frame(x)
  })
  for (k in 1:length(clin.data)) {
    colnames(clin.data[[k]]) <- "label"
  }
  select.group <- c("Elderly", "Young_adult")
  ref.level <- "Young_adult"
  data.type <- "continuous"
  ind.method <- rep("limma", length(data))
  resp.type <- "twoclass"
  paired <- rep(FALSE, length(data))
  meta.res <- MetaDE(data = data, clin.data = clin.data,
    data.type = data.type, resp.type = resp.type,
    response = "label",
    ind.method = ind.method, meta.method = meta.method,
    select.group = select.group, ref.level = ref.level,
    REM.type = "HS",#seed = 2025520,
    paired = paired, tail = tail, parametric = parametric)

  out <- data.frame(ID = allrxns, stringsAsFactors = F,
    MetaAveExpr = rowMeans(flux3.AveExpr, na.rm = T),
    Metalog2FC = rowMeans(flux3.log2FC, na.rm = T),
    MetaPvalue = as.vector(meta.res$meta.analysis$pval),
    MetaFDR = as.vector(meta.res$meta.analysis$FDR),
    MetaZvalue = as.vector(meta.res$meta.analysis$zval))
  out$MetaFDR[is.na(out$MetaFDR)] <- 1
  out$MetaPvalue[is.na(out$MetaPvalue)] <- 1
  out$MetaZvalue[is.na(out$MetaZvalue)] <- 0

  out <- cbind(out, flux3.log2FC)
  out <- cbind(out, flux3.Pvalue)
  out <- cbind(out, flux3.FDR)
  # flux3.annote = Recon3.annote
  # out = cbind(out,flux3.annote[allrxns,])

  # re-calculate pvalue and FDR for pvalue < eps
  idx <- out$MetaPvalue < 1.0
  tmpp <- pnorm(-abs(out$MetaZvalue[idx])) * 2
  tmpfdr <- p.adjust(tmpp, method = "BH")
  out$MetaFDR[idx] <- tmpfdr
  out$MetaPvalue[idx] <- tmpp

  rownames(out) <- allrxns
  tmpfc <- out[, substr(colnames(out), 1, 7) == "log2FC_"]
  colnames(tmpfc) <- gsub("log2FC_", "", colnames(tmpfc))
  out$manyNA <- rowSums(is.na(tmpfc)) >= ncol(tmpfc) / 3
  
  set.seed(2025520)

  return(out)

}

compare_limma_metaDE_eachtissue <- function(DEpro.tissues.lm, DEpro.tissues.MetaDE, outfile, outfile1 = NULL) {

  up1 <- commonMols.filter$FDR005$macaca_uppro_meta
  up2 <- DEpro.tissues.MetaDE$ID[DEpro.tissues.MetaDE$MetaFDR < 0.05 & DEpro.tissues.MetaDE$Metalog2FC > 0.2]
  length(intersect(up1, up2))


  tmpcompare <- list()
  tmpcompare$limma.beta <- DEpro.tissues.lm$MetaLimma[, substr(colnames(DEpro.tissues.lm$MetaLimma), 1, 5) == "beta_"]
  colnames(tmpcompare$limma.beta) <- gsub("beta_", "", colnames(tmpcompare$limma.beta))

  tmpcompare$limma.pval <- DEpro.tissues.lm$MetaLimma[, substr(colnames(DEpro.tissues.lm$MetaLimma), 1, 7) == "Pvalue_"]
  colnames(tmpcompare$limma.pval) <- gsub("Pvalue_", "", colnames(tmpcompare$limma.pval))

  tmpcompare$DE.fc <- DEpro.tissues.MetaDE[, substr(colnames(DEpro.tissues.MetaDE), 1, 7) == "log2FC_"]
  colnames(tmpcompare$DE.fc) <- gsub("beta_", "", colnames(tmpcompare$DE.fc))

  tmpcompare$DE.pval <- DEpro.tissues.MetaDE[, substr(colnames(DEpro.tissues.MetaDE), 1, 7) == "Pvalue_"]
  colnames(tmpcompare$DE.pval) <- gsub("Pvalue_", "", colnames(tmpcompare$DE.pval))

  tmptissues <- colnames(tmpcompare$limma.beta)

  for (i in 1:legnth(tmptissues)) {
    beta1 <- tmpcompare$limma.beta[, 2]
    pval1 <- tmpcompare$limma.pval[, 2]
    names(beta1) <- rownames(tmpcompare$limma.beta)
    names(pval1) <- rownames(tmpcompare$limma.beta)

    fc2 <- tmpcompare$DE.fc[, 2]
    pval2 <- tmpcompare$DE.pval[, 2]
    names(fc2) <- rownames(tmpcompare$DE.fc)
    names(pval2) <- rownames(tmpcompare$DE.fc)
    up1 <- names(beta1)[beta1 > 0 & !is.na(beta1) & pval1 < 0.05]
    up2 <- names(fc2)[fc2 > 0 & !is.na(fc2) & pval2 < 0.05]
    length(up1)
    length(up2)
    length(intersect(up1, up2))



  }




  idx <- intersect(DEpro.tissues.lm$MetaLimma[!DEpro.tissues.lm$MetaLimma$manyNA, ]$ID,
    DEpro.tissues.MetaDE$ID)

  tmpdata <- data.frame(limmabeta = DEpro.tissues.lm$MetaLimma[idx, ]$MetaBeta,
    limmaFDR = DEpro.tissues.lm$MetaLimma[idx, ]$MetaFDR,
    limmaPvalue = DEpro.tissues.lm$MetaLimma[idx, ]$MetaPvalue,
    metalog2FC = DEpro.tissues.MetaDE[idx, ]$Metalog2FC,
    metaFDR = DEpro.tissues.MetaDE[idx, ]$MetaFDR,
    metaPvalue = DEpro.tissues.MetaDE[idx, ]$MetaFDR,
    row.names = idx)
  tmpdata$Sigtype <- (tmpdata$limmaFDR < 0.05) + 0

  pdf(outfile, width = 4.5, height = 4.5)
  tmpcor <- cor.test(tmpdata$limmabeta, tmpdata$metalog2FC)
  rho <- tmpcor$estimate
  pvalue <- tmpcor$p.value

  p <- ggplot(tmpdata, aes(x = limmabeta, y = metalog2FC)) +
    geom_point(aes(color = Sigtype), size = 0.3, alpha = 0.2 + tmpdata$Sigtype * 0.2) +
    geom_density_2d(aes(color = ..level..), size = 0.3, alpha = 0.1) +  # Density contours
    geom_smooth(method = "lm", color = "gray", se = TRUE, alpha = 0.3) +  # Linear regression line
    labs(
      title = "Limma beta vs mean log2FC all",
      x = "Effect size(beta of age)",
      y = "Mean log2FC across 30 tissues(elderly vs young)"
    ) +
    theme_classic() +
    theme(legend.position = "none") + # xlim(c(-3.3,3.3))+ylim(c(-3.3,3.3))+
    annotate("text", x = 0, y = 1, label = paste("Rho =", round(rho, 2)), size = 4, color = "black")
  print(p)

  # only significant
  tmpdata.sig1 <- tmpdata[tmpdata$limmaFDR < 0.05, ]
  tmpcor <- cor.test(tmpdata.sig1$limmabeta, tmpdata.sig1$metalog2FC)
  rho <- tmpcor$estimate
  pvalue <- tmpcor$p.value

  p <- ggplot(tmpdata.sig1, aes(x = limmabeta, y = metalog2FC)) +
    geom_point(aes(color = Sigtype), size = 0.3, alpha = 0.2 + tmpdata.sig1$Sigtype * 0.2) +
    geom_density_2d(aes(color = ..level..), size = 0.3, alpha = 0.1) +  # Density contours
    geom_smooth(method = "lm", color = "gray", se = TRUE, alpha = 0.3) +  # Linear regression line
    labs(
      title = "Limma beta vs mean log2FC significant only",
      x = "Effect size(beta of age)",
      y = "Mean log2FC across 30 tissues(elderly vs young)"
    ) +
    theme_classic() +
    theme(legend.position = "none") + # xlim(c(-3.3,3.3))+ylim(c(-3.3,3.3))+
    annotate("text", x = 0, y = 1, label = paste("Rho =", round(rho, 2)), size = 4, color = "black")
  print(p)
  dev.off()

  # setdiff()
  if (!is.null(outfile1)) {
    # compare beta and log2fc for each tissue
    tissues <- colnames(DEpro.tissues.lm$MetaLimma)[substr(colnames(DEpro.tissues.lm$MetaLimma), 1, 5) == "beta_"]
    tissues <- gsub("beta_", "", tissues)

    idx <- intersect(rownames(DEpro.tissues.lm$MetaLimma), rownames(DEpro.tissues.MetaDE))

    # pdf(file = pdffile,width = 3*5+1,height =3*6)
    # grid.arrange(arrangeGrob(grobs = proteinVoconoPlot,ncol = 5))
    # dev.off()

    pt <- list()
    for (i in 1:length(tissues)) {

      thisname <- tissues[i]
      tmpdata1 <- data.frame(limmabeta = DEpro.tissues.lm$MetaLimma[idx, paste0("beta_", thisname)],
        limmaPvalue = DEpro.tissues.lm$MetaLimma[idx, paste0("Pvalue_", thisname)],
        log2FC = DEpro.tissues.MetaDE[idx, paste0("log2FC_", thisname)],
        dePvalue = DEpro.tissues.MetaDE[idx, paste0("Pvalue_", thisname)],
        row.names = idx)
      tmpdata1$Sigtype <- (tmpdata1$limmaPvalue < 0.05) + 0

      tmpcor <- cor.test(tmpdata1$limmabeta, tmpdata1$log2FC)
      rho <- tmpcor$estimate
      pvalue <- tmpcor$p.value

      pt[[i]] <- ggplot(tmpdata1, aes(x = limmabeta, y = log2FC)) +
        geom_point(aes(color = Sigtype), size = 0.3, alpha = 0.2 + tmpdata1$Sigtype * 0.2) +
        # geom_density_2d(aes(color = ..level..), size = 0.3,alpha = 0.1) +  # Density contours
        # geom_smooth(method = "lm", color = "gray", se = TRUE,alpha = 0.3) +  # Linear regression line
        labs(
          title = thisname,
          x = "Effect size(beta of age)",
          y = "log2FC"
        ) +
        theme_classic() +
        theme(legend.position = "none") + # xlim(c(-3.3,3.3))+ylim(c(-3.3,3.3))+
        annotate("text", x = 0, y = 1, label = paste("Rho =", round(rho, 2)), size = 4, color = "black")

    }

    pdf(file = outfile1, width = 3 * 5 + 1, height = 3 * 6)
    grid.arrange(arrangeGrob(grobs = pt, ncol = 5))
    dev.off()
  }

  # grid.arrange(arrangeGrob(grobs = pt,ncol = 5))

  return(tmpdata)
}


compare_limma_metaDE <- function(DEpro.tissues.lm, DEpro.tissues.MetaDE, 
                                 title, outfile, outfile1 = NULL) {

  idx <- intersect(DEpro.tissues.lm$MetaLimma[!DEpro.tissues.lm$MetaLimma$manyNA, ]$ID,
    DEpro.tissues.MetaDE$ID)

  tmpdata <- data.frame(limmabeta = DEpro.tissues.lm$MetaLimma[idx, ]$MetaBeta,
    limmaFDR = DEpro.tissues.lm$MetaLimma[idx, ]$MetaFDR,
    limmaPvalue = DEpro.tissues.lm$MetaLimma[idx, ]$MetaPvalue,
    metalog2FC = DEpro.tissues.MetaDE[idx, ]$Metalog2FC,
    metaFDR = DEpro.tissues.MetaDE[idx, ]$MetaFDR,
    metaPvalue = DEpro.tissues.MetaDE[idx, ]$MetaFDR,
    row.names = idx)
  tmpdata$Sigtype <- (tmpdata$limmaFDR < 0.05) + 0

  pdf(outfile, width = 3, height = 3)
  tmpcor <- cor.test(tmpdata$limmabeta, tmpdata$metalog2FC)
  rho <- tmpcor$estimate
  pvalue <- tmpcor$p.value

  p <- ggplot(tmpdata, aes(x = limmabeta, y = metalog2FC)) +
    geom_point(aes(color = Sigtype), size = 0.3, alpha = 0.2 + tmpdata$Sigtype * 0.2) +
    geom_density_2d(aes(color = ..level..), size = 0.3, alpha = 0.1) +  # Density contours
    geom_smooth(method = "lm", color = "gray", se = TRUE, alpha = 0.3) +  # Linear regression line
    labs(
      title = paste0(title,  " (all)"),
      x = "Effect size (beta LM)",
      y = "Mean log2FC across 30 tissues"
    ) +
    theme_classic() +
    theme(legend.position = "none") + # xlim(c(-3.3,3.3))+ylim(c(-3.3,3.3))+
    annotate("text", x = 0, y = 1, label = paste("Rho =", round(rho, 3)), size = 4, color = "black")
  print(p)

  # only significant
  tmpdata.sig1 <- tmpdata[tmpdata$limmaFDR < 0.05, ]
  tmpcor <- cor.test(tmpdata.sig1$limmabeta, tmpdata.sig1$metalog2FC)
  rho <- tmpcor$estimate
  pvalue <- tmpcor$p.value

  p <- ggplot(tmpdata.sig1, aes(x = limmabeta, y = metalog2FC)) +
    geom_point(aes(color = Sigtype), size = 0.3, alpha = 0.2 + tmpdata.sig1$Sigtype * 0.2) +
    geom_density_2d(aes(color = ..level..), size = 0.3, alpha = 0.1) +  # Density contours
    geom_smooth(method = "lm", color = "gray", se = TRUE, alpha = 0.3) +  # Linear regression line
    labs(
      title = paste0("Significant ",title),
      x = "Effect size (beta LM)",
      y = "Mean log2FC across 30 tissues"
    ) +
    theme_classic() +
    theme(legend.position = "none") + # xlim(c(-3.3,3.3))+ylim(c(-3.3,3.3))+
    annotate("text", x = 0, y = 1, label = paste("Rho =", round(rho, 3)), size = 4, color = "black")
  print(p)
  
  dev.off()

  # setdiff()
  if (!is.null(outfile1)) {
    # compare beta and log2fc for each tissue
    tissues <- colnames(DEpro.tissues.lm$MetaLimma)[substr(colnames(DEpro.tissues.lm$MetaLimma), 1, 5) == "beta_"]
    tissues <- gsub("beta_", "", tissues)

    idx <- intersect(rownames(DEpro.tissues.lm$MetaLimma), rownames(DEpro.tissues.MetaDE))

    # pdf(file = pdffile,width = 3*5+1,height =3*6)
    # grid.arrange(arrangeGrob(grobs = proteinVoconoPlot,ncol = 5))
    # dev.off()

    pt <- list()
    for (i in 1:length(tissues)) {

      thisname <- tissues[i]
      tmpdata1 <- data.frame(limmabeta = DEpro.tissues.lm$MetaLimma[idx, paste0("beta_", thisname)],
        limmaPvalue = DEpro.tissues.lm$MetaLimma[idx, paste0("Pvalue_", thisname)],
        log2FC = DEpro.tissues.MetaDE[idx, paste0("log2FC_", thisname)],
        dePvalue = DEpro.tissues.MetaDE[idx, paste0("Pvalue_", thisname)],
        row.names = idx)
      tmpdata1$Sigtype <- (tmpdata1$limmaPvalue < 0.05) + 0

      tmpcor <- cor.test(tmpdata1$limmabeta, tmpdata1$log2FC)
      rho <- tmpcor$estimate
      pvalue <- tmpcor$p.value

      pt[[i]] <- ggplot(tmpdata1, aes(x = limmabeta, y = log2FC)) +
        geom_point(aes(color = Sigtype), size = 0.3, alpha = 0.2 + tmpdata1$Sigtype * 0.2) +
        # geom_density_2d(aes(color = ..level..), size = 0.3,alpha = 0.1) +  # Density contours
        # geom_smooth(method = "lm", color = "gray", se = TRUE,alpha = 0.3) +  # Linear regression line
        labs(
          title = thisname,
          x = "Effect size(beta of age)",
          y = "log2FC"
        ) +
        theme_classic() +
        theme(legend.position = "none") + # xlim(c(-3.3,3.3))+ylim(c(-3.3,3.3))+
        annotate("text", x = 0, y = 1, label = paste("Rho =", round(rho, 3)), size = 4, color = "black")

    }

    pdf(file = outfile1, width = 3 * 5 + 1, height = 3 * 6)
    grid.arrange(arrangeGrob(grobs = pt, ncol = 5))
    dev.off()
  }

  # grid.arrange(arrangeGrob(grobs = pt,ncol = 5))

  return(tmpdata)
}


compare_limma_metaDE_ggven <- function(DEpro.tissues.lm, DEpro.tissues.MetaDE,
                                       DEpro.tissues.MetaDE.FEM,
                                 title, outfile) {
  tmplm <- DEpro.tissues.lm$MetaLimma
  tmpxlist  <- list(lm.up = tmplm$ID[tmplm$MetaFDR < 0.05 &
              tmplm$MetaBeta > 0.008 & !tmplm$manyNA],
          lm.down = tmplm$ID[tmplm$MetaFDR < 0.05 &
              tmplm$MetaBeta <  -0.008 & !tmplm$manyNA],
          REM.up = DEpro.tissues.MetaDE$ID[DEpro.tissues.MetaDE$MetaFDR < 0.05 &
              DEpro.tissues.MetaDE$Metalog2FC > 0.2],
          REM.down = DEpro.tissues.MetaDE$ID[DEpro.tissues.MetaDE$MetaFDR < 0.05 &
              DEpro.tissues.MetaDE$Metalog2FC <  -0.2],
          
          FEM.up = DEpro.tissues.MetaDE.FEM$ID[DEpro.tissues.MetaDE.FEM$MetaFDR < 0.05 &
                                                 DEpro.tissues.MetaDE.FEM$Metalog2FC > 0.2],
          FEM.down = DEpro.tissues.MetaDE.FEM$ID[DEpro.tissues.MetaDE.FEM$MetaFDR < 0.05 &
                                                   DEpro.tissues.MetaDE.FEM$Metalog2FC <  -0.2])
  
  pdf(outfile,width = 3.5,height = 3.5)
  
  p <- ggvenn::ggvenn(data = list(LM.up = tmpxlist$lm.up,
                                  Meta.FEM.up = tmpxlist$FEM.up,
                                  Meta.REM.up = tmpxlist$REM.up),
                      fill_color = c("#E64B35FF", "#3C5488FF", "green"),
                      show_percentage = F, stroke_size = 0.5,
                      stroke_alpha = 0.6, text_size = 6)
  print(paste0("UP: (REM & FEM)/REM = ",
        length(intersect(tmpxlist$REM.up,tmpxlist$FEM.up))/length(tmpxlist$REM.up)))
  print(paste0("UP: (REM & LM)/REM = ",
               length(intersect(tmpxlist$REM.up,tmpxlist$lm.up))/length(tmpxlist$REM.up)))
  print(paste0("UP: (FEM & LM)/FEM = ",
               length(intersect(tmpxlist$FEM.up,tmpxlist$lm.up))/length(tmpxlist$FEM.up)))
  
  print(p)
  
  p <- ggvenn::ggvenn(data = list(LM.down = tmpxlist$lm.down,
                                  Meta.FEM.down = tmpxlist$FEM.down,
                                  Meta.REM.down = tmpxlist$REM.down),
                      fill_color = c("#E64B35FF", "#3C5488FF", "green"),
                      show_percentage = F, stroke_size = 0.5,
                      stroke_alpha = 0.6, text_size = 6)
  
  print(paste0("Down: (REM & FEM)/REM = ",
               length(intersect(tmpxlist$REM.down,tmpxlist$FEM.down))/length(tmpxlist$REM.down)))
  print(paste0("Down: (REM & LM)/REM = ",
               length(intersect(tmpxlist$REM.down,tmpxlist$lm.down))/length(tmpxlist$REM.down)))
  print(paste0("Down: (FEM & LM)/FEM = ",
               length(intersect(tmpxlist$FEM.down,tmpxlist$lm.down))/length(tmpxlist$FEM.down)))
  
  print(p)
  
  tmpaa <- list(LM = c(tmpxlist$lm.up,tmpxlist$lm.down),
               Meta.FEM = c(tmpxlist$FEM.up,tmpxlist$FEM.down),
               Meta.REM = c(tmpxlist$REM.up,tmpxlist$REM.down))
  p <- ggvenn::ggvenn(data = tmpaa,
                      fill_color = c("#E64B35FF", "#3C5488FF", "green"),
                      show_percentage = F, stroke_size = 0.5,
                      stroke_alpha = 0.6, text_size = 6)
  
  print(paste0("ALL: (REM & FEM)/REM = ",
               length(intersect(tmpaa$Meta.REM,tmpaa$Meta.FEM))/length(tmpaa$Meta.REM)))
  print(paste0("ALL: (REM & LM)/REM = ",
               length(intersect(tmpaa$Meta.REM,tmpaa$LM))/length(tmpaa$Meta.REM)))
  print(paste0("ALL: (FEM & LM)/FEM = ",
               length(intersect(tmpaa$Meta.FEM,tmpaa$LM))/length(tmpaa$Meta.FEM)))
  
  print(p)
  
  
  dev.off()
  
  # grid.arrange(arrangeGrob(grobs = pt,ncol = 5))
  
  return(p)
}

get_tissue_specific_aging_markers <- function(Metapro, tissue.systems, tissue.color, outfile,
                                              p_cutoff = 0.05, beta_cutoff = 0.008,
                                              topnum = 20,
                                              direction = "UP") {

  probeta <- Metapro[, substr(colnames(Metapro), 1, 5) == "beta_"]
  colnames(probeta) <- gsub("beta_", "", colnames(probeta))

  propvalue <- Metapro[, substr(colnames(Metapro), 1, 7) == "Pvalue_"]
  colnames(propvalue) <- gsub("Pvalue_", "", colnames(propvalue))
  propvalue[is.na(propvalue)] <- 1
  probeta[is.na(probeta)] <- 0
  # probeta[propvalue > 0.05] = 0


  tissue.systems.v <- sort(tissue.systems[names(pro.tissues)])
  ttnames <- names(tissue.systems.v)
  probeta <- probeta[, ttnames]
  propvalue  <- propvalue[, ttnames]

  numsigs <- rowSums(propvalue < 0.05)
  # beta_cutoff
  # top10
  tissue_aging_markers.pro <- list()
  tissue_aging_markers.pro.full <- list()
  maker.tissues <- c()
  for (i in 1:ncol(probeta)) {
    # tmpPvalue = propvalue[,i]
    tmpaa  <- data.frame(ID = rownames(probeta),
      beta = probeta[, i],
      pvalue = propvalue[, i],
      metaFDR = Metapro$MetaFDR,
      numsigs = numsigs,
      row.names = rownames(probeta))
    tmpaa <- tmpaa[order(tmpaa$pvalue), ]

    if (toupper(direction) == "UP") {
      thismarker <- tmpaa$ID[tmpaa$pvalue < p_cutoff & tmpaa$beta  > beta_cutoff &
        tmpaa$numsigs < 3 & tmpaa$metaFDR > 0.05]
    } else {
      thismarker <- tmpaa$ID[tmpaa$pvalue < p_cutoff & tmpaa$beta  < -beta_cutoff &
        tmpaa$numsigs < 3 & tmpaa$metaFDR > 0.05]
    }

    if (length(thismarker) > topnum) {
      thismarker <- thismarker[1:topnum]
    }
    tissue_aging_markers.pro[[i]] <- thismarker
    maker.tissues <- c(maker.tissues, rep(colnames(probeta)[i], length(thismarker)))

  }
  names(tissue_aging_markers.pro) <- colnames(probeta)

  significant_genes <- unlist(tissue_aging_markers.pro) # Combine all significant genes
  expression_matrix <- probeta[significant_genes, , drop = FALSE]

  enbrks <- 0.6 * c(-1, -0.8, -0.4, -0.2, -0.1, -0.05, 0.05, 0.1, 0.2, 0.4, 0.8, 1)
  cname <- colnames(probeta)
  tclass <- data.frame(System = tissue.systems[cname], stringsAsFactors = F, row.names = cname)
  rowclass <- data.frame(System = tissue.systems[maker.tissues], stringsAsFactors = F, row.names = rownames(expression_matrix))
  ann_colors <- list(
    System = tissue.color)

  # Step 1: Identify the top 1 tissue-specific gene for each tissue
  top_genes <- sapply(tissue_aging_markers.pro, function(x) x[1])  # Top gene for each tissue
  # matching_genes <- rownames(expression_matrix) %in% top_genes
  matching_genes <- !duplicated(maker.tissues)
  row_labels <- rep("", nrow(expression_matrix))  # Start with blank labels
  min_vector <- function(x, a) {
    y <- x
    for (i in 1:length(x)) {
      y[i] <- min(x[i], a)
    }
    return(y)
  }
  top_genes_reduce <- substr(top_genes, 1, min_vector(nchar(top_genes), 15))
  row_labels[matching_genes] <- paste0(top_genes_reduce, " (", names(top_genes), ")") # rownames(expression_matrix)[matching_genes]

  # color
  col_fun <- colorRamp2(
    breaks = enbrks,
    colors = colorRampPalette(c("#4DBBD5FF", "gray99", "#E64B35FF"))(length(enbrks))
  )

  # ann col
  col_anno <- HeatmapAnnotation(
    System = tclass$System,
    col = list(System = ann_colors$System),
    annotation_name_side = "left",
    show_annotation_name = TRUE
  )

  # ann row
  row_anno <- rowAnnotation(
    System = rowclass$System,
    col = list(System = ann_colors$System),
    show_annotation_name = FALSE
  )

  # gene markers
  ha_row <- rowAnnotation(
    gene_labels = anno_mark(
      at = which(matching_genes),
      labels = row_labels[matching_genes],
      labels_gp = gpar(fontsize = 10),
      lines_gp = gpar(lwd = 0.5),
      link_width = unit(4, "mm")
    ),
    width = unit(2, "cm") + max_text_width(row_labels[matching_genes])
  )

  # draw
  ht <- Heatmap(
    as.matrix(expression_matrix),
    name = "Beta LM",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_rot = 45,
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 9.5),
    top_annotation = col_anno,
    left_annotation = row_anno,
    right_annotation = ha_row,
    # border_gp = gpar(col = "gray99", lwd = 0.2),
    # rect_gp = gpar(col = "gray99", lwd = 0.2),
    column_title = "Heatmap of aging_beta in 30 tissues",
    heatmap_legend_param = list(
      color_bar = "continuous",
      legend_direction = "vertical"
    )
  )

  # save
  pdf(outfile, width = 12, height = 8)
  draw(ht)
  dev.off()

  out <- list()
  out$makers <- tissue_aging_markers.pro
  out$heatmap <- ht
  return(out)
}

write_tissue_specific_aging_markers <- function(tissue_aging_markers,Metapro,omicstype,excelfile) {
  #omicstype = c('mrna','pro','mets')
  tmpaa <- tissue_aging_markers[[omicstype]]
  markers.up <- tmpaa$up
  markers.down <- tmpaa$down
  # up
  for(i in 1:length(markers.up)){
    tmpname <- names(markers.up)[i]
    thismarkers <- markers.up[[tmpname]]
    tmpcolnames <- c('ID',paste0('beta_',tmpname),paste0('Pvalue_',tmpname))
    tmpfrm <- Metapro[thismarkers,tmpcolnames]
    tmpfrm$Direction <- rep('Increase',nrow(tmpfrm))
    tmpfrm$Tissue <- rep(tmpname,nrow(tmpfrm))
    colnames(tmpfrm) <- c("ID","Beta",'Pvalue','Direction','Tissue')
    tmpfrm[,2:3] <- signif(tmpfrm[,2:3],4)
    if(i == 1){
      markers.up.matrix <- tmpfrm
    }else{
      markers.up.matrix <- rbind(markers.up.matrix,tmpfrm)
    }
  }
  
  #down
  for(i in 1:length(markers.down)){
    tmpname <- names(markers.down)[i]
    thismarkers <- markers.down[[tmpname]]
    tmpcolnames <- c('ID',paste0('beta_',tmpname),paste0('Pvalue_',tmpname))
    tmpfrm <- Metapro[thismarkers,tmpcolnames]
    tmpfrm$Direction <- rep('Decrease',nrow(tmpfrm))
    tmpfrm$Tissue <- rep(tmpname,nrow(tmpfrm))
    colnames(tmpfrm) <- c("ID","Beta",'Pvalue','Direction','Tissue')
    tmpfrm[,2:3] <- signif(tmpfrm[,2:3],4)
    if(i == 1){
      markers.down.matrix <- tmpfrm
    }else{
      markers.down.matrix <- rbind(markers.down.matrix,tmpfrm)
    }
  }
  
  markers.all.matrix <- list()
  markers.all.matrix[["Up_markers"]] <- markers.up.matrix
  markers.all.matrix[["Down_markers"]] <- markers.down.matrix
  openxlsx::write.xlsx(markers.all.matrix, file = excelfile, keepNA = T)
  
  return(markers.all.matrix)
}

get_Metacorrlation_met2pro <- function(met.tissues, pro.tissues.v, Metapro, inputname) {
  # Hypoxanthine
  teron <- list()
  teron[[length(met.tissues)]] <- ""
  # inputname = "Hypoxanthine"
  xnames <- c()
  for (i in 1:length(met.tissues)) {
    tmp <- met.tissues[[i]]
    if (any(rownames(tmp) == inputname)) {
      teron[[i]] <- met.tissues[[i]][inputname, ]
      xnames <- c(xnames, names(met.tissues)[i])
    }
  }
  names(teron) <- names(met.tissues)


  # get_correlation protein
  corTeron2pro <- get_corr(teron, pro.tissues.v)
  corTeron2pro.Pvalue <- str2num(list_element_select(corTeron2pro, rownames(Metapro), "Pvalue"))
  corTeron2pro.cor <- str2num(list_element_select(corTeron2pro, rownames(Metapro), "cor"))
  corTeron2pro.Pvalue[is.na(corTeron2pro.Pvalue)] <- 1

  x <- list(p = corTeron2pro.Pvalue)
  xx.pro <- MetaDE.pvalue(x = x, meta.method = "AW")
  tmpdata <- data.frame(ID = rownames(corTeron2pro.Pvalue),
    log2FC = rowMeans(corTeron2pro.cor, na.rm = T),
    Pvalue = as.vector(xx.pro$meta.analysis$FDR))
  # tmpdata = del_duplicate_rows_forpro(tmpdata,pro.whole.nofilter.header)
  # tmpdata$ID = tmpdata$gene
  rownames(tmpdata) <- tmpdata$ID
  tmpdata$Pvalue[tmpdata$Pvalue < 1e-30] <- 1e-30
  return(tmpdata)

}


volcono_plot_species <- function(Metamrna, title, FCcutoff = 0.005) {
  tmpdata <- Metamrna
  tmpdata$Pvalue <- tmpdata$MetaFDR
  tmpdata$log2FC <- tmpdata$MetaBeta
  tmpdata$log2FC[tmpdata$Pvalue == 1] <- 0

  idx <- tmpdata$manyNA == FALSE
  nstudy <- sum(substr(colnames(tmpdata), 1, 7) == "Pvalue_")
  p <- plot_DEflux(tmpdata[idx, ], alpha = 0.2, num.showlab = 5, FCcutoff = FCcutoff,
    xlab = paste0("LM beta (across ", nstudy, " tissues)", collapse = ""),
    fixpointsize = 1, labSize = 6,
    ylab = bquote(~ -Log[10] ~ italic(FDR)),
    title = title) + #+ xlim(c(-1,1))
    theme(axis.line = element_line(size = 1),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16)
    )
  return(p)
}

get_slid_mrnas_3_species <- function(MetageneMouse, Metamrna, MetageneGETX) {

  out <- list()
  betacut.mouse <- abs(sort(MetageneMouse$MetaBeta)[nrow(MetageneMouse) * 0.1]) # 0.1/18months
  betacut.macaca <- abs(sort(Metamrna$MetaBeta)[nrow(Metamrna) * 0.1]) # 0.1/24years
  betacut.human <- abs(sort(MetageneGETX$MetaBeta)[nrow(MetageneGETX) * 0.1]) # 0.1/60years

  # betacut.mouse = 0 #0.2/21months
  # betacut.macaca = 0 #0.2/27years
  # betacut.human = 0 #0.2/80years

  slids <- (2:12) * 100
  for (i in slids) {
    tmpaa  <- MetageneMouse[order(MetageneMouse$MetaPvalue), ]
    mouse_upmrna <- tmpaa$ID[tmpaa$MetaBeta > betacut.mouse & tmpaa$manyNA == FALSE][1:i]
    mouse_downmrna <- tmpaa$ID[tmpaa$MetaBeta < -betacut.mouse & tmpaa$manyNA == FALSE][1:i]


    tmpaa  <- Metamrna[order(Metamrna$MetaPvalue), ]
    macaca_upmrna <- tmpaa$ID[tmpaa$MetaBeta > betacut.macaca & tmpaa$manyNA == FALSE][1:i]
    macaca_downmrna <- tmpaa$ID[tmpaa$MetaBeta < -betacut.macaca & tmpaa$manyNA == FALSE][1:i]


    tmpaa  <- MetageneGETX[order(MetageneGETX$MetaPvalue), ]
    human_upmrna <- tmpaa$ID[tmpaa$MetaBeta > betacut.human & tmpaa$manyNA == FALSE][1:i]
    human_downmrna <- tmpaa$ID[tmpaa$MetaBeta < -betacut.human & tmpaa$manyNA == FALSE][1:i]
    out[[paste0("TOP", i)]] <- list(mouse_upmrna = mouse_upmrna,
      mouse_downmrna = mouse_downmrna,
      macaca_upmrna = macaca_upmrna,
      macaca_downmrna = macaca_downmrna,
      human_upmrna = human_upmrna,
      human_downmrna = human_downmrna)
  }
  return(out)

}

compare_GO_C2 <- function(commonMols.filter, TERM2GENE, outfile) {

  out.rho.up <- rep(0, 3)
  names(out.rho.up)  <- c("C2_macaca_vs_mouse", "C2_macaca_vs_human", "C2_human_vs_mouse")
  out.rho.down <- out.rho.up

  mouse_upmrna  <- GOenrichment.C2(commonMols.filter$FDR005$mouse_upmrna, TERM2GENE)
  mouse_downmrna  <- GOenrichment.C2(commonMols.filter$FDR005$mouse_downmrna, TERM2GENE)
  macaca_upmrna  <- GOenrichment.C2(commonMols.filter$FDR005$macaca_upmrna_meta, TERM2GENE)
  macaca_downmrna  <- GOenrichment.C2(commonMols.filter$FDR005$macaca_downmrna_meta, TERM2GENE)
  human_upmrna  <- GOenrichment.C2(commonMols.filter$FDR005$human_upmrna, TERM2GENE)
  human_downmrna  <- GOenrichment.C2(commonMols.filter$FDR005$human_downmrna, TERM2GENE)
  tmpGO <- list(mouse_upmrna = mouse_upmrna,
    mouse_downmrna = mouse_downmrna,
    macaca_upmrna = macaca_upmrna,
    macaca_downmrna = macaca_downmrna,
    human_upmrna = human_upmrna,
    human_downmrna = human_downmrna)
  gotype <- c("C2")
  tt <- 0
  for (j in 1:length(gotype)) {

    k1 <- tmpGO$mouse_upmrna[tmpGO$mouse_upmrna$type == gotype[j], ]
    k2 <- tmpGO$macaca_upmrna[tmpGO$macaca_upmrna$type == gotype[j], ]
    k3 <- tmpGO$human_upmrna[tmpGO$human_upmrna$type == gotype[j], ]

    tt <- tt + 1
    vid <- intersect(rownames(k1), rownames(k2))
    out.rho.up[tt] <- cor.test(log2(k1[vid, ]$ratio), log2(k2[vid, ]$ratio))$estimate

    tt <- tt + 1
    vid <- intersect(rownames(k2), rownames(k3))
    out.rho.up[tt] <- cor.test(log2(k2[vid, ]$ratio), log2(k3[vid, ]$ratio))$estimate

    tt <- tt + 1
    vid <- intersect(rownames(k1), rownames(k3))
    out.rho.up[tt] <- cor.test(log2(k1[vid, ]$ratio), log2(k3[vid, ]$ratio))$estimate
  }

  ## go down
  # rho up
  tt <- 0
  for (j in 1:length(gotype)) {

    k1 <- tmpGO$mouse_downmrna[tmpGO$mouse_downmrna$type == gotype[j], ]
    k2 <- tmpGO$macaca_downmrna[tmpGO$macaca_downmrna$type == gotype[j], ]
    k3 <- tmpGO$human_downmrna[tmpGO$human_downmrna$type == gotype[j], ]

    tt <- tt + 1
    vid <- intersect(rownames(k1), rownames(k2))
    out.rho.down[tt] <- cor.test(log2(k1[vid, ]$ratio), log2(k2[vid, ]$ratio))$estimate

    tt <- tt + 1
    vid <- intersect(rownames(k2), rownames(k3))
    out.rho.down[tt] <- cor.test(log2(k2[vid, ]$ratio), log2(k3[vid, ]$ratio))$estimate

    tt <- tt + 1
    vid <- intersect(rownames(k1), rownames(k3))
    out.rho.down[tt] <- cor.test(log2(k1[vid, ]$ratio), log2(k3[vid, ]$ratio))$estimate
  }
  out.rho <- list(up = out.rho.up, down = out.rho.down)
  out <- list(rho = out.rho,
    GO = tmpGO)

  # plot
  C2_compare <- out
  vids <- intersect(rownames(C2_compare$GO$mouse_upmrna), rownames(C2_compare$GO$macaca_upmrna))
  vids <- intersect(vids, rownames(C2_compare$GO$human_upmrna))
  combined_data_ratio.up <- data.frame(
    pathway = vids,
    log2ratio_mouse = log2(C2_compare$GO$mouse_upmrna[vids, ]$ratio),
    log2ratio_human = log2(C2_compare$GO$human_upmrna[vids, ]$ratio),
    log2ratio_macaca = log2(C2_compare$GO$macaca_upmrna[vids, ]$ratio),
    Sigtype = (C2_compare$GO$human_upmrna[vids, ]$p.adjust < 0.05) + 0,
    Sigtypemcc = (C2_compare$GO$macaca_upmrna[vids, ]$p.adjust < 0.05) + 0
  )

  vids <- intersect(rownames(C2_compare$GO$mouse_downmrna), rownames(C2_compare$GO$macaca_downmrna))
  vids <- intersect(vids, rownames(C2_compare$GO$human_downmrna))
  combined_data_ratio.down <- data.frame(
    pathway = vids,
    log2ratio_mouse = log2(C2_compare$GO$mouse_downmrna[vids, ]$ratio),
    log2ratio_human = log2(C2_compare$GO$human_downmrna[vids, ]$ratio),
    log2ratio_macaca = log2(C2_compare$GO$macaca_downmrna[vids, ]$ratio),
    Sigtype = (C2_compare$GO$human_downmrna[vids, ]$p.adjust < 0.05) + 0,
    Sigtypemcc = (C2_compare$GO$macaca_downmrna[vids, ]$p.adjust < 0.05) + 0
  )

  tmpfun_plot_c2 <- function(combined_data_ratio.up, dir = "up") {
    p <- list()
    # human vs macaca
    rho <- cor.test(combined_data_ratio.up$log2ratio_macaca, combined_data_ratio.up$log2ratio_human, method = "spearman")$estimate
    p[[1]] <- ggplot(combined_data_ratio.up, aes(x = log2ratio_human, y = log2ratio_macaca)) +
      geom_point(size = 0.3, alpha = 0.4) +
      geom_density_2d(aes(color = ..level..), size = 0.3, alpha = 0.2) +  # Density contours
      # geom_smooth(method = "lm", color = "black", se = FALSE) +  # Linear regression line
      labs(
        title = paste0("Human_vs_macaque ", dir, "-regulated"),
        x = "log2ratio (Human)",
        y = "log2ratio (Macaque)"
      ) +
      theme_bw() +
      theme(legend.position = "none") + # xlim(c(-3.3,3.3))+ylim(c(-3.3,3.3))+
      annotate("text", x = 0, y = 5, label = paste("Rho =", round(rho, 2)), size = 4, color = "black")

    # human vs mouse
    rho <- cor.test(combined_data_ratio.up$log2ratio_mouse, combined_data_ratio.up$log2ratio_human, method = "spearman")$estimate
    p[[2]] <- ggplot(combined_data_ratio.up, aes(x = log2ratio_human, y = log2ratio_mouse)) +
      geom_point(size = 0.3, alpha = 0.4) +
      geom_density_2d(aes(color = ..level..), size = 0.3, alpha = 0.2) +  # Density contours
      # geom_smooth(method = "lm", color = "black", se = FALSE) +  # Linear regression line
      labs(
        title = paste0("Human_vs_mouse ", dir, "-regulated"),
        x = "log2ratio (Human)",
        y = "log2ratio (Mouse)"
      ) +
      theme_bw() +
      theme(legend.position = "none") + # xlim(c(-3.3,3.3))+ylim(c(-3.3,3.3))+
      annotate("text", x = 0, y = 5, label = paste("Rho =", round(rho, 2)), size = 4, color = "black")

    # macaca vs mouse
    rho <- cor.test(combined_data_ratio.up$log2ratio_mouse, combined_data_ratio.up$log2ratio_macaca, method = "spearman")$estimate
    p[[3]] <- ggplot(combined_data_ratio.up, aes(x = log2ratio_macaca, y = log2ratio_mouse)) +
      geom_point(size = 0.3, alpha = 0.4) +
      geom_density_2d(aes(color = ..level..), size = 0.3, alpha = 0.2) +  # Density contours
      # geom_smooth(method = "lm", color = "black", se = FALSE) +  # Linear regression line
      labs(
        title = paste0("Macaque_vs_mouse ", dir, "-regulated"),
        x = "log2ratio (Macaque)",
        y = "log2ratio (Mouse)"
      ) +
      theme_bw() +
      theme(legend.position = "none") + # xlim(c(-3.3,3.3))+ylim(c(-3.3,3.3))+
      annotate("text", x = 0, y = 5, label = paste("Rho =", round(rho, 2)), size = 4, color = "black")
    return(p)

  }

  pdf(outfile, width = 3.5, height = 3.5)
  p1 <- tmpfun_plot_c2(combined_data_ratio.up)
  print(p1[[1]])
  print(p1[[2]])
  print(p1[[3]])

  p2 <- tmpfun_plot_c2(combined_data_ratio.down, "down")
  print(p2[[1]])
  print(p2[[2]])
  print(p2[[3]])
  dev.off()
  out$plots.up <- p1
  out$plots.down <- p2

  return(out)
}



compare_GO <- function(commonMols.filter) {

  out.rho.up <- rep(0, 6)
  names(out.rho.up)  <- c("KEGG_macaca_vs_mouse", "KEGG_macaca_vs_human", "KEGG_human_vs_mouse",
    "BP_macaca_vs_mouse", "BP_macaca_vs_human", "BP_human_vs_mouse")
  out.rho.down <- out.rho.up

  mouse_upmrna  <- GOenrichment(commonMols.filter$FDR005$mouse_upmrna)
  mouse_downmrna  <- GOenrichment(commonMols.filter$FDR005$mouse_downmrna)
  macaca_upmrna  <- GOenrichment(commonMols.filter$FDR005$macaca_upmrna_meta)
  macaca_downmrna  <- GOenrichment(commonMols.filter$FDR005$macaca_downmrna_meta)
  human_upmrna  <- GOenrichment(commonMols.filter$FDR005$human_upmrna)
  human_downmrna  <- GOenrichment(commonMols.filter$FDR005$human_downmrna)
  tmpGO <- list(mouse_upmrna = mouse_upmrna,
    mouse_downmrna = mouse_downmrna,
    macaca_upmrna = macaca_upmrna,
    macaca_downmrna = macaca_downmrna,
    human_upmrna = human_upmrna,
    human_downmrna = human_downmrna)
  gotype <- c("KEGG", "BP")
  tt <- 0
  for (j in 1:length(gotype)) {

    k1 <- tmpGO$mouse_upmrna[tmpGO$mouse_upmrna$type == gotype[j], ]
    k2 <- tmpGO$macaca_upmrna[tmpGO$macaca_upmrna$type == gotype[j], ]
    k3 <- tmpGO$human_upmrna[tmpGO$human_upmrna$type == gotype[j], ]

    tt <- tt + 1
    vid <- intersect(rownames(k1), rownames(k2))
    out.rho.up[tt] <- cor.test(log2(k1[vid, ]$ratio), log2(k2[vid, ]$ratio))$estimate

    tt <- tt + 1
    vid <- intersect(rownames(k2), rownames(k3))
    out.rho.up[tt] <- cor.test(log2(k2[vid, ]$ratio), log2(k3[vid, ]$ratio))$estimate

    tt <- tt + 1
    vid <- intersect(rownames(k1), rownames(k3))
    out.rho.up[tt] <- cor.test(log2(k1[vid, ]$ratio), log2(k3[vid, ]$ratio))$estimate
  }

  ## go down
  # rho up
  tt <- 0
  for (j in 1:length(gotype)) {

    k1 <- tmpGO$mouse_downmrna[tmpGO$mouse_downmrna$type == gotype[j], ]
    k2 <- tmpGO$macaca_downmrna[tmpGO$macaca_downmrna$type == gotype[j], ]
    k3 <- tmpGO$human_downmrna[tmpGO$human_downmrna$type == gotype[j], ]

    tt <- tt + 1
    vid <- intersect(rownames(k1), rownames(k2))
    out.rho.down[tt] <- cor.test(log2(k1[vid, ]$ratio), log2(k2[vid, ]$ratio))$estimate

    tt <- tt + 1
    vid <- intersect(rownames(k2), rownames(k3))
    out.rho.down[tt] <- cor.test(log2(k2[vid, ]$ratio), log2(k3[vid, ]$ratio))$estimate

    tt <- tt + 1
    vid <- intersect(rownames(k1), rownames(k3))
    out.rho.down[tt] <- cor.test(log2(k1[vid, ]$ratio), log2(k3[vid, ]$ratio))$estimate
  }
  out.rho <- list(up = out.rho.up, down = out.rho.down)
  return(out.rho)


}

get_GO_slid_mrnas <- function(slid_mrnas) {
  slids <- (2:12) * 100

  out <- list()

  for (i in slids) {
    xid <- paste0("TOP", i)
    print(xid)
    tmpmrna <- slid_mrnas[[xid]]
    mouse_upmrna  <- GOenrichment(tmpmrna$mouse_upmrna)
    mouse_downmrna  <- GOenrichment(tmpmrna$mouse_downmrna)
    macaca_upmrna  <- GOenrichment(tmpmrna$macaca_upmrna)
    macaca_downmrna  <- GOenrichment(tmpmrna$macaca_downmrna)
    human_upmrna  <- GOenrichment(tmpmrna$human_upmrna)
    human_downmrna  <- GOenrichment(tmpmrna$human_downmrna)
    out[[xid]] <- list(mouse_upmrna = mouse_upmrna,
      mouse_downmrna = mouse_downmrna,
      macaca_upmrna = macaca_upmrna,
      macaca_downmrna = macaca_downmrna,
      human_upmrna = human_upmrna,
      human_downmrna = human_downmrna)
  }

  return(out)
}


plot_GOenrich_ratio <- function(slid_mrnas, slid_GO, outpath) {
  #### overlap genes
  slids <- (2:12) * 100
  k <- 0
  out.rho.up <- matrix(0, 6, length(slids))
  colnames(out.rho.up)  <-  paste0("TOP", slids)
  rownames(out.rho.up)  <- c("KEGG_macaca_vs_mouse", "KEGG_macaca_vs_human", "KEGG_human_vs_mouse",
    "BP_macaca_vs_mouse", "BP_macaca_vs_human", "BP_human_vs_mouse")

  out.rho.down <- out.rho.up
  for (i in slids) {
    xid <- paste0("TOP", i)
    tmpmrna <- slid_mrnas[[xid]]
    tmpGO <- slid_GO[[xid]]

    p <- ggvenn::ggvenn(data = list(mrna.up.human = tmpmrna$human_upmrna,
      mrna.up.macaca = tmpmrna$macaca_upmrna,
      mrna.up.mouse = tmpmrna$mouse_upmrna),
    fill_color = c("#E64B35FF", "#3C5488FF", "green"),
    show_percentage = F, stroke_size = 0.5,
    stroke_alpha = 0.6, text_size = 9)
    pdf(paste0(outpath, "/", xid, "_ggvenn_mrna_up.pdf"), width = 8)
    print(p)
    dev.off()

    p <- ggvenn::ggvenn(data = list(mrna.down.human = tmpmrna$human_downmrna,
      mrna.down.macaca = tmpmrna$macaca_downmrna,
      mrna.down.mouse = tmpmrna$mouse_downmrna),
    fill_color = c("#E64B35FF", "#3C5488FF", "green"),
    show_percentage = F, stroke_size = 0.5,
    stroke_alpha = 0.6, text_size = 9)
    pdf(paste0(outpath, "/", xid, "_ggvenn_mrna_down.pdf"), width = 8)
    print(p)
    dev.off()

    k <- k + 1
    # rho up
    gotype <- c("KEGG", "BP")
    tt <- 0
    for (j in 1:length(gotype)) {

      k1 <- tmpGO$mouse_upmrna[tmpGO$mouse_upmrna$type == gotype[j], ]
      k2 <- tmpGO$macaca_upmrna[tmpGO$macaca_upmrna$type == gotype[j], ]
      k3 <- tmpGO$human_upmrna[tmpGO$human_upmrna$type == gotype[j], ]

      tt <- tt + 1
      vid <- intersect(rownames(k1), rownames(k2))
      out.rho.up[tt, k] <- cor.test(log2(k1[vid, ]$ratio), log2(k2[vid, ]$ratio))$estimate

      tt <- tt + 1
      vid <- intersect(rownames(k2), rownames(k3))
      out.rho.up[tt, k] <- cor.test(log2(k2[vid, ]$ratio), log2(k3[vid, ]$ratio))$estimate

      tt <- tt + 1
      vid <- intersect(rownames(k1), rownames(k3))
      out.rho.up[tt, k] <- cor.test(log2(k1[vid, ]$ratio), log2(k3[vid, ]$ratio))$estimate
    }

    ## go down
    # rho up
    tt <- 0
    for (j in 1:length(gotype)) {

      k1 <- tmpGO$mouse_downmrna[tmpGO$mouse_downmrna$type == gotype[j], ]
      k2 <- tmpGO$macaca_downmrna[tmpGO$macaca_downmrna$type == gotype[j], ]
      k3 <- tmpGO$human_downmrna[tmpGO$human_downmrna$type == gotype[j], ]

      tt <- tt + 1
      vid <- intersect(rownames(k1), rownames(k2))
      out.rho.down[tt, k] <- cor.test(log2(k1[vid, ]$ratio), log2(k2[vid, ]$ratio))$estimate

      tt <- tt + 1
      vid <- intersect(rownames(k2), rownames(k3))
      out.rho.down[tt, k] <- cor.test(log2(k2[vid, ]$ratio), log2(k3[vid, ]$ratio))$estimate

      tt <- tt + 1
      vid <- intersect(rownames(k1), rownames(k3))
      out.rho.down[tt, k] <- cor.test(log2(k1[vid, ]$ratio), log2(k3[vid, ]$ratio))$estimate
    }

  }
  out.rho <- list(up = out.rho.up, down = out.rho.down)
  return(out.rho)
}


plot_heatmap_fromMetadata <- function(Metamrna, threecomm, filename = filename, tissue.systems = tissue.systems,
                                      height = 6, width = 8, cluster_cols = F, cluster_rows = T,
                                      show_rownames = T,
                                      show_colnames = T,
                                      fontsize_row = 11, fontsize_col = 11) {
  # mcc mrna
  mrnafc.common <- Metamrna[threecomm,
    substr(colnames(Metamrna), 1, 7) == "log2FC_"]
  colnames(mrnafc.common) <- gsub("log2FC_", "", colnames(mrnafc.common))

  mrnapval.common <- Metamrna[threecomm,
    substr(colnames(Metamrna), 1, 7) == "Pvalue_"]
  colnames(mrnapval.common) <- gsub("Pvalue_", "", colnames(mrnapval.common))


  mrnapval.common <- mrnapval.common[rownames(mrnafc.common), colnames(mrnafc.common)]

  #
  display_matrix <- matrix(" ", nrow(mrnapval.common), ncol(mrnapval.common))
  display_matrix[mrnapval.common < 0.1] <- "."
  display_matrix[mrnapval.common >= 0.01 & mrnapval.common < 0.05] <- "*"
  display_matrix[mrnapval.common >= 0.001 & mrnapval.common < 0.01] <- "**"
  display_matrix[mrnapval.common < 0.001] <- "***"

  enbrks <- c(-1, -0.8, -0.4, -0.2, -0.05, -0.02, 0.02, 0.05, 0.2, 0.4, 0.8, 1)
  # DAscore.all_1 = DAscore.all
  # DAscore.all_1[DAscore.all_1 == 0] = NA

  cname <- colnames(mrnafc.common)

  tclass <- data.frame(System = tissue.systems[cname], row.names = cname)
  mrnafc.common.v <- mrnafc.common
  mrnafc.common.v[mrnafc.common.v > 1] <- 1
  mrnafc.common.v[mrnafc.common.v < -1] <- -1

  ann_colors <- list(
    System = tissue.color)

  heatmap_mrna <- pheatmap::pheatmap(mrnafc.common.v, scale = "none", cluster_rows = cluster_rows,
    cluster_cols = cluster_cols, annotation_colors = ann_colors,
    fontsize_row = fontsize_row, fontsize_col = fontsize_col,
    annotation_col = tclass, show_rownames = show_rownames,
    show_colnames = show_rownames,
    breaks = enbrks,
    angle_col = 45,
    treeheight_row = 20, treeheight_col = 20, legend = T,
    display_numbers = display_matrix,
    # color=colorRampPalette(c('#3B4992','gray95','#BB0021'))(9),
    # colorRampPalette(c('#008280','gray95','#BB0021'))(11),
    color = colorRampPalette(c("#4DBBD5FF", "gray95", "#E64B35FF"))(11),
    # color=colorRampPalette(c('#3B4992','gray99','#EE0000'))(11),
    # color=colorRampPalette(c('#008280FF','gray99','#E64B35FF'))(11),
    file = filename,
    height = height, width = width)
  return(heatmap_mrna)
}



runclustering <- function(d1, title = title, maxK = 6, ngene = dim(d1)[1], distance = "pearson", clusterAlg = "km") {
  # distance can be "euclidean"
  require(ConsensusClusterPlus)
  mads <- apply(d1, 1, mad)
  d1 <- d1[rev(order(mads))[1:ngene], ]
  d1 <- sweep(d1, 1, apply(d1, 1, median, na.rm = T))
  if (dir.exists(title)) {
    a <- 1
  } else {
    dir.create(title)
  }
  # title<-tempdir()
  b <- as.matrix(d1)
  results <- ConsensusClusterPlus(b, maxK = maxK, reps = 1000, pItem = 0.8, pFeature = 1, title = title, clusterAlg = clusterAlg, distance = distance, seed = 2025520, plot = "pdf")
  return(results)
}


get_trj_data <- function(promet.tissues, promet.tissues.info,
                         promet.tissues.Z, mfuzz.promet.whole,
                         tissue.systems, tissue.color, outfile, center = 1) {
  tissue_trajectory <- data.frame()
  tissue_names <- names(promet.tissues)
  for (k in 1:8) {
    tmpGeneList <- names(mfuzz.promet.whole$cluster)[mfuzz.promet.whole$cluster == k]
    for (i in 1:length(promet.tissues)) {

      M1 <- promet.tissues[[i]]
      idgene1 <- intersect(tmpGeneList, rownames(M1))
      M1 <- M1[idgene1, ]
      metadata.tissue <- promet.tissues.info[[i]]
      M1.median <- t(aggregate(t(M1), by = list(metadata.tissue$stage),
        FUN = median, na.rm = T))
      M1.median.Z <- standardise_matrix(M1.median)
      mean_x <- colMeans(M1.median.Z, na.rm = T)
      mean_x <- mean_x - mean_x[center]

      tmpdata2 <- data.frame(expr = mean_x, stringsAsFactors = F,
        stage = c(1, 2, 3, 4),
        group = factor(c("Juvenile", "Young", "Middle_aged", "Elderly"),
          level = c("Juvenile", "Young", "Middle_aged", "Elderly")),
        cluster = rep(k, 4),
        tissue = rep(tissue_names[i], 4)
      )
      if (nrow(tissue_trajectory) < 1) {
        tissue_trajectory <- tmpdata2
      } else {
        tissue_trajectory <- rbind(tissue_trajectory, tmpdata2)
      }
    }
  }

  tissue_trajectory$tissue_system <- tissue.systems[tissue_trajectory$tissue]
  tissue_trajectory$tissue_system_color <- tissue.color[tissue_trajectory$tissue_system]
  # to matrix
  tissue_trajectory_matrix <- matrix(0, nrow(tissue_trajectory) / length(tissue_names), length(tissue_names))
  for (i in 1:length(tissue_names)) {
    tmptr <- tissue_trajectory[tissue_trajectory$tissue == tissue_names[i], ]
    tissue_trajectory_matrix[, i] <- tmptr$expr
  }
  rownames(tissue_trajectory_matrix) <- paste(tmptr$group, tmptr$cluster, sep = "_C")
  colnames(tissue_trajectory_matrix) <- tissue_names

  ## MAA
  tissues <- names(promet.tissues.Z)
  clusterdist <- matrix(0, length(tissues), 8)
  clusteramplitude <- matrix(0, length(tissues), 8)
  clusteramplitude_xx <- matrix(0, length(tissues), 8)
  for (i in 1:length(tissues)) {
    mstd <-   promet.tissues.Z[[i]]
    for (j in 1:8) {
      tgene <- names(mfuzz.promet.whole$cluster)[mfuzz.promet.whole$cluster == j]
      tgene <- intersect(tgene, rownames(mstd))
      bx <- t(t(mstd[tgene, ]) - mfuzz.promet.whole$centers[j, ])
      clusterdist[i, j] <- mean(sqrt(rowSums(bx^2) / (ncol(bx) - 1)), na.rm = T)
      if (length(tgene) < 2) {
        clusteramplitude[i, j] <- NA
        next
      }
      tamp <- colMeans(mstd[tgene, ], na.rm = T)
      clusteramplitude[i, j] <- abs(tamp[4] - tamp[center])
      clusteramplitude_xx[i, j] <- tamp[4] - tamp[center]
    }
  }
  rownames(clusteramplitude) <- tissues
  rownames(clusteramplitude_xx) <- tissues
  rownames(clusterdist) <- tissues


  ## plot figure 5B tissue trajactory
  tissue_tr_plot <- list()
  for (i in 1:8) {
    idx <- tissue_trajectory$cluster == i
    tissue_tr_plot [[i]] <- ggplot(tissue_trajectory[idx, ],
      aes(x = group, y = expr, group  = tissue)) +
      geom_line(aes(color = tissue_system), alpha = 0.6,
        position = position_dodge(0.2), size = 2) + scale_color_npg() +
      lghplot.addtheme(hjust = 1, size = 14) + ggtitle(paste0("Cluster ", i)) +
      theme(axis.line = element_line(size = 1.3)) + xlab("") + ylab("") +
      scale_y_continuous(labels = scales::comma_format(accuracy = 0.1))
  }
  pdf(file = outfile, height = 7, width = 11)
  grid.arrange(arrangeGrob(grobs = tissue_tr_plot, ncol = 4, heights = c(4, 4),
    top = textGrob("Average molecular aging trajectory in 30 tissues",
      gp = gpar(fontface = "bold",  fontsize = 20)),
    left = textGrob("Z values", gp = gpar(fontface = "bold",  fontsize = 18), rot = 90)))
  dev.off()


  out <- list()
  out$tissue_trajectory <- tissue_trajectory
  out$tissue_trajectory_matrix <- tissue_trajectory_matrix
  out$MAA <- clusteramplitude_xx
  out$plot <- tissue_tr_plot
  return(out)
}
