
# A multi-omics molecular landscape of 30 tissues in aging female rhesus macaques (Macaca mulatta)

**Code author**: Gong-Hua Li 
  
**Project authors**: Gong-Hua Li#, Xiang-Qing Zhu#, Fu-Hui Xiao#, Xilong Zhao, Longbao Lv, Fan-Qian Yin, Le Chang, Ming-Xia Ge, Qiang Wang, Jing Zhao, Chuan Tian, Zian Li, Guangping Ruan, Rongqing Pang, Jing Gao, Lihua Ma, Xing-Hua Pan*, Qing-Peng Kong*  
  
**Project contact**: Qing-Peng Kong (kongqp@mail.kiz.ac.cn) or Xing-Hua Pan (xinghuapan@aliyun.com) 
  
**Code & data analysis contact**: Gong-Hua Li (ligonghua@mail.kiz.ac.cn)  
  
**Last Update**: 2025-03-20    

---

## 1. Data Availability  
The omics data are freely and publicly available. Level-2 data, including transcriptome raw counts and log2 peak areas for the proteome and metabolome, are accessible on Figshare1.  Level-1 raw data, including FASTQ files for the transcriptome and raw LC-MS/MS data for the proteome and metabolome, have been deposited in the OMIX, China National Center for Bioinformation / Beijing Institute of Genomics, Chinese Academy of Sciences (https://ngdc.cncb.ac.cn/omix: accession no. OMIX001777 for transcriptome, OMIX001778 for proteome, and OMIX001779 for metabolome). 

Li, G.-H. & Kong, Q.-P. A multi-omics molecular landscape of 30 tissues in aging female rhesus macaques (Macaca mulatta). figshare Dataset, https://doi.org/10.6084/m6089.figshare.26963386 (2024).
---

## 2. Data Description  
*(Content to be added)*  

---

## 3. Install Required Packages

```r
# R Packages from CRAN
packages_cran <- c("ggplot2", "ggbreak", "reshape2", "Hmisc", "pracma", "car", "ggrepel", 
                   "grid", "gridExtra", "rlist", "ggsci", "scales", "data.table", "openxlsx",
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
