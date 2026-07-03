## RNA-seq project packages
## Run the installation section only when setting up a new R library/server.
## Run the loading section every R session.

## 1) Install packages if needed

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

cran_pkgs <- c(
  "tidyverse",
  "pheatmap",
  "RColorBrewer",
  "patchwork",
  "writexl",
  "msigdbr"
)

bioc_pkgs <- c(
  "DESeq2",
  "edgeR",
  "limma",
  "apeglm",
  "EnhancedVolcano",
  "BiocParallel",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "fgsea",
  "SummarizedExperiment"
)

to_install_cran <- cran_pkgs[
  !vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(to_install_cran) > 0) {
  install.packages(to_install_cran, repos = "https://cloud.r-project.org")
}

to_install_bioc <- bioc_pkgs[
  !vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(to_install_bioc) > 0) {
  BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)
}


## 2) Load packages for RNA-seq analysis

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(apeglm)
  library(EnhancedVolcano)
  library(BiocParallel)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(fgsea)
  library(msigdbr)
  library(SummarizedExperiment)
  library(pheatmap)
  library(RColorBrewer)
  library(patchwork)
  library(writexl)
  library(grid)
})
