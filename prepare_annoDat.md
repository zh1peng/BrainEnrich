---
title: "How to prepare annoDat for the analysis"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{prepare_annoDat}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




``` r
# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")}

  bioc_packages <- c("org.Hs.eg.db", "clusterProfiler", "DOSE","ReactomePA", "reactome.db")
  for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg)
    }
  }

build_Anno <- getFromNamespace("build_Anno", "DOSE")
```
# Gene ontology (GO) gene sets

``` r
get_GO_data <- getFromNamespace("get_GO_data", "clusterProfiler")
annoData <- get_GO_data("org.Hs.eg.db", ont = 'BP', keytype = "SYMBOL")
#saveRDS(annoData, file = "inst/extdata/geneSets/GO_BP.rds")
annoData <- get_GO_data("org.Hs.eg.db", ont = 'MF', keytype = "SYMBOL")
#saveRDS(annoData, file = "inst/extdata/geneSets/GO_MF.rds")
annoData <- get_GO_data("org.Hs.eg.db", ont = 'CC', keytype = "SYMBOL")
save(annoData, file = "inst/extdata/geneSets/GO_CC.rds")
# check version of org.Hs.eg.db
packageVersion("clusterProfiler")
packageVersion("org.Hs.eg.db")
```

# DisGeNET gene sets

``` r
 get_DGN_data <- getFromNamespace("get_DGN_data", "DOSE")
 annoData <- get_DGN_data()
 save(annoData, file = "inst/extdata/geneSets/DGN.rds")
 packageVersion("DOSE")
```

# KEGG

``` r
prepare_KEGG <- getFromNamespace("prepare_KEGG", "clusterProfiler")
annoData <- prepare_KEGG("hsa", "MKEGG", "kegg")
#saveRDS(annoData, file = "inst/extdata/geneSets/KEGG.rds")
packageVersion("clusterProfiler")
```

# WikiPathways

``` r
prepare_WP_data <- getFromNamespace("prepare_WP_data", "clusterProfiler")
wpdata <- prepare_WP_data("Homo sapiens")
TERM2GENE <- wpdata$WPID2GENE
TERM2NAME <- wpdata$WPID2NAME
build_Anno <- getFromNamespace("build_Anno", "DOSE")
annoData <- build_Anno(TERM2GENE, TERM2NAME)
#saveRDS(annoData, file = "inst/extdata/geneSets/WikiPathways.rds")
packageVersion("clusterProfiler")
```

# Reactome

``` r
get_Reactome_DATA <- getFromNamespace("get_Reactome_DATA", "ReactomePA")
annoData <- get_Reactome_DATA("human")
#saveRDS(annoData, file = "inst/extdata/geneSets/Reactome.rds")
```

# SynGO (how to make annoDat from xlsx)

``` r
# Define URL for downloading SynGO data
url <- "https://www.syngoportal.org/data/SynGO_bulk_download_release_20231201.zip"
# Create temporary file and directory
zip_path <- tempfile(fileext = ".zip")
temp_dir <- tempdir()
# Clean-up function to ensure temp files are removed
on.exit({
  if (file.exists(zip_path)) unlink(zip_path)
  if (file.exists(extracted_file)) unlink(extracted_file)
}, add = TRUE)
extracted_file <- file.path(temp_dir, "syngo_ontologies.xlsx")
download.file(url, zip_path, mode = "wb")
unzip(zip_path, files = "syngo_ontologies.xlsx", exdir = temp_dir)


# Read the data from the extracted file
data <- read_xlsx(extracted_file)

# Process TERM2GENE
TERM2GENE <- data %>%
  dplyr::select(id, hgnc_symbol) %>%
  dplyr::rename(term = id, gene = hgnc_symbol) %>%
  dplyr::mutate(gene = strsplit(gene, ", ")) %>%
  unnest(cols = c(gene))

# Process TERM2NAME
TERM2NAME <- data %>%
  dplyr::select(id, name) %>%
  dplyr::rename(term = id, description = name)

# Use build_Anno function from DOSE package
build_Anno <- getFromNamespace("build_Anno", "DOSE")
annoData <- build_Anno(TERM2GENE, TERM2NAME)
#saveRDS(annoData, file = "inst/extdata/geneSets/SynGO.rds")
```

# CellTypes Seidlitz2020 (how to prepare annoData with a csv file)
### Seidlitz, J., Nadig, A., Liu, S., Bethlehem, R. A., Vértes, P. E., Morgan, S. E., ... & Raznahan, A. (2020). Transcriptomic and cellular decoding of regional brain vulnerability to neurogenetic disorders. Nature communications, 11(1), 3358.

``` r
# Define URL for Seidlitz2020
url <- "https://github.com/jms290/PolySyn_MSNs/blob/master/Data/AHBA/celltypes_PSP.csv?raw=true"

# Create temporary file
temp_file <- tempfile()
# Clean-up function to ensure temp files are removed
on.exit({
  if (file.exists(temp_file)) unlink(temp_file)
}, add = TRUE)

download.file(url, temp_file, mode = "wb")

# Read CSV file
TERM2GENE <- read.csv(temp_file) %>%
  dplyr::mutate(term = class) %>%
  filter(gene != "") %>%
  dplyr::select(term, gene)
TERM2NAME <- TERM2GENE %>%
  dplyr::mutate(description = term) %>%
  dplyr::select(term, description)
annoData <- build_Anno(TERM2GENE, TERM2NAME)
#saveRDS(annoData, file = "inst/extdata/geneSets/CellTypes_Seidlitz2020.rds")
```


# CellTypes Lake2018 (how to prepare annoData with a gmt file)
### Lake, B. B., Chen, S., Sos, B. C., Fan, J., Kaeser, G. E., Yung, Y. C., ... & Zhang, K. (2018). Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain. Nature biotechnology, 36(1), 70-80.

``` r
# Define URL for Lake2018
url <- "https://github.com/molecular-neuroimaging/Imaging_Transcriptomics/raw/main/imaging_transcriptomics/data/geneset_LAKE.gmt"
# Create temporary file
temp_file <- tempfile()
# Clean-up function to ensure temp files are removed
on.exit({
  if (file.exists(temp_file)) unlink(temp_file)
}, add = TRUE)


download.file(url, temp_file, mode = "wb")
# Read GMT file
TERM2GENE <- suppressWarnings({
  clusterProfiler::read.gmt(temp_file) %>%
    filter(gene != "") %>%
    dplyr::select(term, gene)
})
TERM2NAME <- TERM2GENE %>%
  dplyr::mutate(description = term) %>%
  dplyr::select(term, description)
annoData <- build_Anno(TERM2GENE, TERM2NAME)
#saveRDS(annoData, file = "inst/extdata/geneSets/CellTypes_Lake2018.rds")
```

# CellTypes Martins2021 (how to prepare annoData with a gmt file)
### Imaging transcriptomics: Convergent cellular, transcriptomic, and molecular neuroimaging signatures in the healthy adult human brain. Daniel Martins, Alessio Giacomel, Steven CR Williams, Federico Turkheimer, Ottavia Dipasquale, Mattia Veronese, PET templates working group. Cell Reports.

``` r
# Define URL for Martins2021
url <- "https://github.com/molecular-neuroimaging/Imaging_Transcriptomics/raw/main/imaging_transcriptomics/data/geneset_Pooled.gmt"
# Create temporary file
temp_file <- tempfile()
# Clean-up function to ensure temp files are removed
on.exit({
  if (file.exists(temp_file)) unlink(temp_file)
}, add = TRUE)
download.file(url, temp_file, mode = "wb")
# Read GMT file
TERM2GENE <- suppressWarnings({
  clusterProfiler::read.gmt(temp_file) %>%
    filter(gene != "") %>%
    dplyr::select(term, gene)
})
TERM2NAME <- TERM2GENE %>%
  dplyr::mutate(description = term) %>%
  dplyr::select(term, description)
annoData <- build_Anno(TERM2GENE, TERM2NAME)
#saveRDS(annoData, file = "inst/extdata/geneSets/CellTypes_Martins2021.rds")
```



# render me

