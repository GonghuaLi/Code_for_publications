
# A multi-omics molecular landscape of 30 tissues in aging female rhesus macaques (Macaca mulatta)

**Code author**: Gong-Hua Li 
  
**Project authors**: Gong-Hua Li#, Xiang-Qing Zhu#, Fu-Hui Xiao#, Xilong Zhao, Longbao Lv, Fan-Qian Yin, Le Chang, Ming-Xia Ge, Qiang Wang, Jing Zhao, Chuan Tian, Zian Li, Guangping Ruan, Rongqing Pang, Jing Gao, Lihua Ma, Xing-Hua Pan*, Qing-Peng Kong*  
  
**Project contact**: Qing-Peng Kong (kongqp@mail.kiz.ac.cn) or Xing-Hua Pan (xinghuapan@aliyun.com) 
  
**Code & data analysis contact**: Gong-Hua Li (ligonghua@mail.kiz.ac.cn)  
  
**Last Update**: 2025-03-20    

---

## 1. Data Availability  
The omics data are freely and publicly available. 

Level-2 data, including transcriptome raw counts and log2 peak areas for the proteome and metabolome, are accessible on figshare (https://doi.org/10.6084/m6089.figshare.26963386).  

Level-1 raw data, including FASTQ files for the transcriptome and raw LC-MS/MS data for the proteome and metabolome, have been deposited in the OMIX, China National Center for Bioinformation / Beijing Institute of Genomics, Chinese Academy of Sciences (https://ngdc.cncb.ac.cn/omix: accession no. OMIX001777 for transcriptome, OMIX001778 for proteome, and OMIX001779 for metabolome). 

---

## 2. Data Description  
The following table provides an overview description of the datasets used in this project:

| Data_name                                      | Data_type       | Description                                              |
|------------------------------------------------|-----------------|----------------------------------------------------------|
| rnaData.Rdata                                  | R transcriptome | Macaque 30 tissue raw counts data                        |
| pro.tissues_solid_tissues.Rdata                | R proteome      | Macaque 30 tissue proteome data (log2 peak area)         |
| pro.whole.fdr0.01_from_NOVO_remap_solid_tissues.Rdata | R proteome      | Macaque whole body proteome data (log2 peak area)        |
| met.header.all.hmdb_curated.Rdata              | R metabolome    | Macaque metabolite information data                      |
| met.tissues_solid_tissues.Rdata                | R metabolome    | Macaque 30 tissue metabolome data (log2 peak area)       |
| met_whole_from_novo_solid_tissues.Rdata        | R metabolome    | Macaque whole body metabolome data (log2 peak area)      |
| Macaca_mulatta.Mmul_10.99.ensemble_symbol_biotype.txt | txt gene_information | Macaque gene information file                     |
| GETX.cpm_21_Female_tissues.RData               | R transcriptome | Human female 21 tissues log2 TMM-normalized CPM data from GTEX |
| geneInfo_encodev38.tab                         | txt gene_information | Human gene information file                       |
| Mus_agingNature_reads_list.Rdata               | R transcriptome | Mouse 17 tissue raw counts data from Schaum et al., Nature, 2020 |
| mouse_geneInfo_GRCm39v113.txt                  | txt gene_information | Mouse gene information file                       |
| mouse_proteome_aging.xlsx                      | xlsx proteome   | Mouse aging-related proteins from Keele et al., Cell reports, 2023 |

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

### Change the rootdir for the code file "Code_for_Macaca_30tissue_aging.Rmd"
please change the line 27:

**rootdir <- "/home/ligh/pubproject/MCMT/"**   to  your own path, such as **rootdir <- "/your/path/"**

**Note:** in your path, a "data" subdirectory is required for reproduce our results, which contains all data that is available at figshare ((https://doi.org/10.6084/m6089.figshare.26963386).


### Generate HTML Report
```bash
Rscript -e "rmarkdown::render('Code_for_Macaca_30tissue_aging.Rmd', output_file='Code_to_Results.html')"
```


### Generate PDF Report
```bash
Rscript -e "rmarkdown::render('Code_for_Macaca_30tissue_aging.Rmd', output_format='pdf_document', output_file='Code_to_Results.pdf')"
```
