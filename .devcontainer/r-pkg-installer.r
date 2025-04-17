options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install(
    version = "devel",
    ask     = FALSE,
    update  = FALSE
)

install.packages(
    c(
        "tidyverse",
        "tidymodels",
        "kableExtra",
        "languageserver"
    ),
    dependencies = TRUE
)

BiocManager::install(
    c(
        "GEOquery",
        "affy",
        "limma",
        "hgu133plus2.db",
        "tximport",
        "edgeR",
        "rhdf5",
        "ALDEx2"
    ),
    ask = FALSE,
    update = FALSE
)
