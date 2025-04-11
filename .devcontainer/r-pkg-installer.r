if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install(version = "3.21")
BiocManager::install("GEOquery")
BiocManager::install("affy")
BiocManager::install("limma")
BiocManager::install("hgu133plus2.db")
BiocManager::install("tximport")
BiocManager::install("edgeR")

install.packages("tidyverse")
install.packages("tidymodels")
install.packages("kableExtra")
install.packages("languagesever")
