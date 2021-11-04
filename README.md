# Strongyloides RNA-seq Browser Preprocessing and Analysis 
Preprocessing and analysis related to the [*Strongyloides* RNA-seq Browser](https://hallemlab.shinyapps.io/Strongyloides_RNAseq_Browser/), a web-based Shiny App for browsing and on-demand analysis of *Strongyloides spp.* RNA-seq datasets.

## Table of Contents  
1. [General Information](#general-information)
2. [App Access](#app-access)
3. [Sources](#sources)
4. [License](#license)
5. [Authors](#authors)
6. [Update Notes](#update-notes)

## General Information
This repository contains non-responsive code for the pre-processing and analysis for *Strongyloides spp* RNA-seq datasets. Preprocessed data is used as inputs for the *Strongyloides* RNA-seq Browser, a Shiny Web App for on-demand browsing and analysis of published bulk *Strongyloides* RNA-seq data.  It also contains analysis code used to generate results discussed in Bryant, DeMarco, and Hallem (2021). [A preprint version of the manuscript is available on bioRxiv.](https://www.biorxiv.org/content/10.1101/2021.02.18.431867v1). A final version of this manuscript is [published in G3](https://pubmed.ncbi.nlm.nih.gov/33823530/).

The sections below describe the contents of the primary subfolders within this repository.

### Preprocessing  
This folder contains RMarkdown files for each species; these files contain code for the alignment of raw reads, data filtering and normalization, voom variance-normalization of count data, and collection of gene annotation information. RMarkdown files generate outputs that act as essential inputs to the *Strongyloides* RNA-seq Browser App. RMarkdown files are knitted into PDF files; see those PDF files for plots illustrating the effect of filtering and normallization on raw data inputs. See the readme file located in the subfolder for additional details. This folder also contains data files used in offline data pre-processing, including study design files and raw transcripts per million datasets for each species, as well as an Ensembl Compara database of parasite gene sets.  

### Analysis  
This folder contains RMarkdown files for each species, in which we present example analyses. These analyses include hierarchical clustering and principal component analyses of species samples, limma-voom differential expression (example results echoing analyses completed by the Shiny App), benchmarking of Browser datasets and differential expression analyses against previously published datasets, and an example gene set enrichment analysis. Rmd files are knitted into html files, and both PCA and Benchmarking plots are saved in a Plot subfolder in the Outputs folder. See the readme file located in the subfolder for additional details.  

### Supplemental Materials  
This folder contains supplemental figures and files from the manuscript of Bryant, DeMarco, and Hallem (2021). See subfolder README for details.  

### Time_course_DGE  
This folder contains an RMarkdown file and cache that tests for differential expression on longitudinal data using the ImpulseDE2 package. The ImpulseDE2 method contrasts with DE algorithms such as limma that treat time points independently. Instead, ImpulseDE2 seeks to identify genes that display specific trajectories of differential gene expression over time. An Impulse model is designed to capture 4 different expression trajectories: monotonous decrease, monotonous increase, transient decrease (valley), and transient increase (peak). The analysis included in this folder specifically tests an implementation of ImpulseDE2 analysis on a subset of the *S. stercoralis* dataset: life stages FLF, PF, iL3, iL3a.  

## App Access
To access a stable deployment of the *Strongyloides* RNA-seq Browser Web App, please visit: [hallemlab.shinyapps.io/strongyloides_rnaseq_browser/](https://hallemlab.shinyapps.io/Strongyloides_RNAseq_Browser/)  

To view full source code for the *Strongyloides* RNA-seq Browser, please visit the [app repository](https://github.com/HallemLab/Strongyloides_RNAseq_Browser). 

## Sources
* [Shiny](https://shiny.rstudio.com/) - UI framework
* [kallisto](https://pachterlab.github.io/kallisto/) - Lightweight RNA-seq pseudoalignment method
* [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) - Diferential gene expression
* [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) - Gene set enrichment analysis
* *Strongyloides* RNA-seq datasets:
  - [Stolzfus *et al* 2012](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0001854)
  - [Akimori *et al* 2021](https://www.nature.com/articles/s41598-021-87478-3)
  - [Hunt *et al* 2016](https://www.nature.com/articles/ng.3495)
  - [Hunt *et al* 2018](https://www.nature.com/articles/s41598-018-23514-z)
* [WormbaseParasite](https://parasite.wormbase.org/index.html) - Gene annotations and reference transcriptomes
* [DIYTranscriptomics](http://diytranscriptomics.com/) - Virtual asynchronous course where the authors learned best practices for RNA-seq data analysis; provided primary pipeline for data pre-processing and analysis

## License  
This project is licensed under the MIT License. 

## Authors  
* [Astra Bryant, PhD](https://github.com/astrasb)
* [Stephanie DeMarco, PhD](https://github.com/sfdemarco)
* [Elissa Hallem, PhD](https://github.com/ehallem)

## Update Notes
Updated by A.S.B. on approximately 11-3-21 to include *S. stercoralis* free-living male RNA-seq data published by Akimori *et al* 2021. 