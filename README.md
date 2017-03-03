# Workflow link

Link to complete workflow http://bioinformaticsfmrp.github.io/TCGAWorkflow/

# Abstract

Biotechnological advances in sequencing have led to an explosion of
    publicly available data via large international consortia such as [The
    Cancer Genome Atlas (TCGA)](http://cancergenome.nih.gov/), [The
    Encyclopedia of DNA Elements (ENCODE)](http://www.encodeproject.org/),
    and [The NIH Roadmap Epigenomics Mapping Consortium
    (Roadmap)](http://www.roadmapepigenomics.org/). These projects have
    provided unprecedented opportunities to interrogate the epigenome of
    cultured cancer cell lines as well as normal and tumor tissues with high
    genomic resolution. The [Bioconductor](http://www.bioconductor.org/)
    project offers more than 1,000 open-source software and statistical
    packages to analyze high-throughput genomic data. However, most packages
    are designed for specific data types (e.g. expression, epigenetics,
    genomics) and there is no one comprehensive tool that provides a
    complete integrative analysis of the resources and data provided by all
    three public projects. A need to create an integration of these
    different analyses was recently proposed. In this workflow, we provide a
    series of biologically focused integrative analyses of different
    molecular data. We describe how to download, process and prepare TCGA
    data and by harnessing several key Bioconductor packages, we describe
    how to extract biologically meaningful genomic and epigenomic data.
    Using Roadmap and ENCODE data, we provide a work plan to identify
    biologically relevant functional epigenomic elements associated with
    cancer. 
    
To illustrate our workflow, we analyzed two types of brain
tumors: low-grade glioma (LGG) versus high-grade glioma (glioblastoma
multiform or GBM). 

All the  package landing pages used in this workflow can be found through the [biocViews interface](http://www.bioconductor.org/packages/release/BiocViews.html#___Software).
    
**Keywords:** Epigenomics, Genomics, Cancer, non-coding, TCGA, ENCODE, Roadmap, Bioinformatics.

# Installation

To be able to execute all the steps of this workflow please install it with the following code:

```R
source("https://bioconductor.org/biocLite.R")
deps <- c("pathview","clusterProfiler","ELMER", "DO.db","GO.db", 
          "ComplexHeatmap","EDASeq", "TCGAbiolinks","AnnotationHub",
          "gaia","ChIPseeker","minet","BSgenome.Hsapiens.UCSC.hg19",
          "MotifDb","MotIV", "rGADEM", "motifStack","RTCGAToolbox")
for(pkg in deps)  if (!pkg %in% installed.packages()) biocLite(pkg, dependencies = TRUE)
deps <- c("devtools","DT","pbapply","readr","circlize")
for(pkg in deps)  if (!pkg %in% installed.packages())  install.packages(pkg,dependencies = TRUE)
devtools::install_github("BioinformaticsFMRP/TCGAWorkflowData")
devtools::install_github("BioinformaticsFMRP/TCGAWorkflow", dependencies = TRUE)
```

# Docker image

A complete enviroment with all packages installed is available as Docker image, which can be easily run on Mac OS, Windows and Linux systems. The image can be obtained from Docker Hub: https://hub.docker.com/r/tiagochst/tcgabiolinksgui/

Download image:
```{bash, eval = FALSE}
docker pull tiagochst/tcgabiolinksgui
```

To run R from the command line:
```{bash, eval = FALSE}
docker run -ti tiagochst/tcgabiolinksgui R
```

To run RStudio Server (user: rstudio, password: rstudio): 
```{bash, eval = FALSE}
docker run -p 8787:8787 tiagochst/tcgabiolinksgui
```

For more information please check: https://docs.docker.com/

# Loading packages

At the beginning of each section, the packages required to execute the code will be loaded. However the following packages are required for all sections.

- TCGAWorkflowData: this package contains the data necessary to execute each of the analysis
steps. This is a subset of the downloaded to make the example faster. For a real analysis, please
use all the data available.
- DT: we will use it to visualize the results

```R
library(TCGAWorkflowData)
library(DT)
```

# Building the vignette
```R
devtools::build_vignettes()
```
