#setwd('/home/ligh/projects/X/projects/mecfs/')
setwd('E:/iHuman/projects/mecfs')
library(rGPMM)
library(data.table)


plot_DEflux1 <- function (DEflux, pcutoff = 0.05, num.showlab = 10, title = "DEflux", 
          FCcutoff = 0, alpha = 1, xlab = bquote(~Log[2] ~ "fold change"), addinglab = NULL,
          fixpointsize = NULL, ylab = bquote(~-Log[10] ~ italic(P)), 
          labSize = 4, legendLabSize = 12) 
{
  require(EnhancedVolcano)
  res2 = DEflux
  keyvals <- rep("gray50", nrow(res2))
  names(keyvals) <- rep("NS", nrow(res2))
  keyvals[which(res2$log2FC > FCcutoff & res2$Pvalue < pcutoff)] <- "Brown"
  names(keyvals)[which(res2$log2FC > FCcutoff & res2$Pvalue < 
                         pcutoff)] <- "UP"
  keyvals[which(res2$log2FC < -FCcutoff & res2$Pvalue < pcutoff)] <- "darkblue"
  names(keyvals)[which(res2$log2FC < -FCcutoff & res2$Pvalue < 
                         pcutoff)] <- "Down"
  if (!is.element("ID", colnames(res2))) {
    res2$ID = rownames(res2)
  }
  else {
    rownames(res2) = res2$ID
  }
  thelab = res2$ID
  uptmp = res2[res2$log2FC > 0, ]
  downtmp = res2[res2$log2FC < 0, ]
  xid = sort.int(uptmp$Pvalue, decreasing = F, index.return = T)$ix
  uptmp = uptmp[xid, ]
  if (nrow(uptmp) > num.showlab) {
    uplab = uptmp$ID[1:num.showlab]
  }
  else {
    uplab = uptmp$ID
  }
  xid = sort.int(downtmp$Pvalue, decreasing = F, index.return = T)$ix
  downtmp = downtmp[xid, ]
  if (nrow(uptmp) > num.showlab) {
    downlab = downtmp$ID[1:num.showlab]
  }
  else {
    downlab = downtmp$ID
  }
  ymax = max(-log10(res2$Pvalue)) + 1
  thelab = c(uplab, downlab)
  vindx = res2[thelab, ]$Pvalue < pcutoff
  thelab = thelab[vindx]
  thelab = unique(c(thelab,addinglab))
  all_label = res2$ID
  all_label[!is.element(all_label, thelab)] = NA
  if (is.null(fixpointsize)) {
    p = EnhancedVolcano(res2, ylim = c(0, ymax), lab = all_label, 
                        x = "log2FC", y = "Pvalue", title = title, border = "full", 
                        titleLabSize = 18, FCcutoff = FCcutoff, cutoffLineWidth = 0, 
                        cutoffLineType = "blank", axisLabSize = 18, subtitle = NULL, 
                        cutoffLineCol = "white", gridlines.minor = F, gridlines.major = F, 
                        xlab = xlab, ylab = ylab, pCutoff = pcutoff, colCustom = keyvals, 
                        colAlpha = 4/5, legendPosition = "bottom", legendLabSize = legendLabSize, 
                        legendIconSize = 3, drawConnectors = TRUE, widthConnectors = 0.5, 
                        pointSize = -alpha * log10(res2$Pvalue), labSize = labSize, 
                        colConnectors = "black", max.overlaps = Inf)
  }
  else {
    p = EnhancedVolcano(res2, ylim = c(0, ymax), lab = all_label, 
                        x = "log2FC", y = "Pvalue", title = title, border = "full", 
                        titleLabSize = 18, FCcutoff = FCcutoff, cutoffLineWidth = 0, 
                        cutoffLineType = "blank", axisLabSize = 18, subtitle = NULL, 
                        cutoffLineCol = "white", gridlines.minor = F, gridlines.major = F, 
                        xlab = xlab, ylab = ylab, pCutoff = pcutoff, colCustom = keyvals, 
                        colAlpha = 4/5, legendPosition = "bottom", legendLabSize = legendLabSize, 
                        legendIconSize = 3, drawConnectors = TRUE, widthConnectors = 0.5, 
                        pointSize = fixpointsize, labSize = labSize, colConnectors = "black", 
                        max.overlaps = Inf)
  }
  return(p)
}


#DEflux
deflux = file2frame('./results/GSE245661/DEflux.txt')
rownames(deflux) = deflux$ID
idx = deflux$subSystemes == "Alanine and aspartate metabolism" & deflux$Pvalue < 0.05
addinglab = deflux$ID[idx]
addinglab = c(addinglab, '4HGLSDm','PHCHGSm', 
              'G6PDH2r', 'PGL', 'RPE', 'RPI', 'TKT1', 'TKT2', 'GND')
xa = deflux$ID
xa[!is.element(deflux$ID,addinglab)] = NA

pdf('./results/GSE245661/DEflux_large_size1f.pdf')
plot_DEflux1(DEflux = deflux,num.showlab = 2, alpha = 1.5,#addinglab = addinglab,
             FCcutoff = 0.2,labSize = 6)+
  lghplot.addtheme(size = 20,sizex = 18,sizey = 18)+ 
  ggrepel::geom_label_repel(aes(label = xa),size = 4,force = 20,box.padding = 0.1,max.overlaps = 1000,force_pull = 1)
dev.off()

WDA.bootstrap = file2frame('./results/GSE245661/WDA.cen.bootstrap.txt')
idaa = WDA.bootstrap$subsystem != "Exchange/demand reaction" & 
  WDA.bootstrap$fdr<0.05 & abs(WDA.bootstrap$DAscore) > 0.2 
pdf("./results/GSE245661/DAscore_remove_nonsig.pdf")
plot_DAscore(WDA.bootstrap$DAscore[idaa], WDA.bootstrap$subsystem[idaa], 
             enlarge = 1, gtitle = "Enzymatic and transport reactions", 
             WDA.bootstrap$fdr[idaa], gsize = 14, gsizex = 12, gsizey = 12)
dev.off()





# key metabolites
load('./results/GSE245661/keymetabolites.Rdata')
