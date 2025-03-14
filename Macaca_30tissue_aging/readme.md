
# A multi-omics molecular landscape of 30 tissues in aging female rhesus macaques (Macaca mulatta)

**Code author**: Gong-Hua Li 
  
**Project authors**: Gong-Hua Li#, Xiang-Qing Zhu#, Fu-Hui Xiao#, Xilong Zhao, Longbao Lv, Fan-Qian Yin, Le Chang, Ming-Xia Ge, Qiang Wang, Jing Zhao, Chuan Tian, Zian Li, Guangping Ruan, Rongqing Pang, Jing Gao, Lihua Ma, Xing-Hua Pan*, Qing-Peng Kong*  
  
**Project contact**: Qing-Peng Kong (kongqp@mail.kiz.ac.cn) or Xing-Hua Pan (xinghuapan@aliyun.com) 
  
**Code & data analysis contact**: Gong-Hua Li (ligonghua@mail.kiz.ac.cn)  
  
**Last Update**: 2025-03-20    

---

## 1. Data Availability  
*(Content to be added)*  

---

## 2. Data Description  
*(Content to be added)*  

---

## 3. Install Required Packages

```r
# R Packages from CRAN
packages_cran <- c("ggplot2", "ggbreak", "reshape2", "Hmisc", "pracma", "car", "ggrepel", 
                   "grid", "gridExtra", "rlist", "ggsci", "scales", "data.table", 
                   "RColorBrewer", "readxl", "circlize", "ggvenn")

for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}


# Bioconductor Packages
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages_bioc <- c("M3C", "sva", "edgeR", "limma", "Mfuzz", "EnhancedVolcano", 
                   "pheatmap", "org.Hs.eg.db", "clusterProfiler", "fgsea", "msigdbr", 
                   "enrichplot", "ComplexHeatmap", "ConsensusClusterPlus", "samr", "DESeq2")

for (pkg in packages_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}


# Install MetaDE (GitHub)
if (!require("MetaDE", quietly = TRUE)) {
  install.packages("remotes")
  remotes::install_github("metaOmics/MetaDE")
}
```
---

## 4. Run the RMarkdown Analysis

### Generate HTML Report
```bash
Rscript -e "rmarkdown::render('Code_for_Macaca_30tissue_aging.Rmd', output_file='Code_to_Results.html')"
```

### Generate PDF Report
```bash
Rscript -e "rmarkdown::render('Code_for_Macaca_30tissue_aging.Rmd', output_format='pdf_document', output_file='Code_to_Results.pdf')"
```
