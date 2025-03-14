
# A multi-omics molecular landscape of 30 tissues in aging female rhesus macaques (Macaca mulatta)

**Author**: Gong-Hua Li  
**Last Update**: 2024-06-20  

---

## 1. Install Required Packages

### R Packages from CRAN
```r
packages_cran <- c("ggplot2", "ggbreak", "reshape2", "Hmisc", "pracma", "car", "ggrepel", 
                   "grid", "gridExtra", "rlist", "ggsci", "scales", "data.table", 
                   "RColorBrewer", "readxl", "circlize", "ggvenn")

for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
```

### Bioconductor Packages
```r
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
```

### Install MetaDE (GitHub)
```r
if (!require("MetaDE", quietly = TRUE)) {
  install.packages("remotes")
  remotes::install_github("metaOmics/MetaDE")
}
```

---

## 2. Data Availability  
*(Content to be added)*  

---

## 3. Data Description  
*(Content to be added)*  

---

## 4. Run the RMarkdown Analysis

### Generate HTML Report
```bash
Rscript -e "rmarkdown::render('Code_for_Macaca_30tissue_aging.Rmd')"
```

### Generate PDF Report
```bash
Rscript -e "rmarkdown::render('Code_for_Macaca_30tissue_aging.Rmd', output_format='pdf_document')"
```
