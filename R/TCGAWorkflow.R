#'  Biotechnological advances in sequencing have led to an explosion of
#'  publicly available data via large international consortia such as The
#'  Cancer Genome Atlas (TCGA), The
#'  Encyclopedia of DNA Elements (ENCODE),
#'  and The NIH Roadmap Epigenomics Mapping Consortium
#'  (Roadmap). These projects have provided unprecedented opportunities to interrogate the epigenome of
#'  cultured cancer cell lines as well as normal and tumor tissues with high
#'  genomic resolution. The Bioconductor project offers more than 1,000 open-source software and statistical
#'  packages to analyze high-throughput genomic data. However, most packages
#'  are designed for specific data types (e.g. expression, epigenetics,
#'  genomics) and there is no one comprehensive tool that provides a
#'  complete integrative analysis of the resources and data provided by all
#'  three public projects. A need to create an integration of these
#'  different analyses was recently proposed. In this workflow, we provide a
#'  series of biologically focused integrative analyses of different
#'  molecular data. We describe how to download, process and prepare TCGA
#'  data and by harnessing several key Bioconductor packages, we describe
#'  how to extract biologically meaningful genomic and epigenomic data.
#'  Using Roadmap and ENCODE data, we provide a work plan to identify
#'  biologically relevant functional epigenomic elements associated with
#'  cancer. 
#'  To illustrate our workflow, we analyzed two types of brain tumors: 
#'  low-grade glioma (LGG) versus high-grade glioma (glioblastoma multiform or GBM). 
#' @docType package
#' @name TCGAWorkflow
#' @import ELMER  downloader SummarizedExperiment  TCGAWorkflowData knitr
#' @import gaia ChIPseeker ComplexHeatmap
#' @import minet c3net biomaRt pathview
#' @import  motifStack TCGAbiolinks graphics
#' @import circlize pbapply GenomeInfoDb ggplot2 ggthemes parallel
#' @importFrom rGADEM GADEM nOccurrences getPWM consensus nMotifs
#' @importFrom MotIV viewAlignments motifMatch 
#' @importFrom AnnotationHub query AnnotationHub
#' @importFrom clusterProfiler bitr
#' @import DT BSgenome.Hsapiens.UCSC.hg19 GenomicRanges RTCGAToolbox
NULL
