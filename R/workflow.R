#' Biotechnological advances in sequencing have led to an explosion of publicly 
#' available data via large international consortia such as The Cancer Genome Atlas (TCGA), 
#' The Encyclopedia of DNA Elements (ENCODE), and The NIH Roadmap Epigenomics Mapping Consortium (Roadmap).
#'  These projects have provided unprecedented opportunities to interrogate the 
#'  epigenome of cultured cancer cell lines as well as normal and tumor tissues 
#'  with high genomic resolution. The Bioconductor project offers more than 1,000 
#'  open-source software and statistical packages to analyze high-throughput genomic data. 
#'  However, most packages are designed for specific data types (e.g. expression, epigenetics, genomics) 
#'  and there is no one comprehensive tool that provides a complete integrative analysis 
#'  of the resources and data provided by all three public projects. 
#'  A need to create an integration of these different analyses was recently proposed. 
#'  In this workflow, we provide a series of biologically focused integrative 
#'  analyses of different molecular data. We describe how to download, 
#'  process and prepare TCGA data and by harnessing several key Bioconductor packages, 
#'  we describe how to extract biologically meaningful genomic and epigenomic data. 
#'  Using Roadmap and ENCODE data, we provide a work plan to identify 
#'  biologically relevant functional epigenomic elements associated with cancer.
#'  To illustrate our workflow, we analyzed two types of brain tumors : 
#'  low-grade glioma (LGG) versus high-grade glioma (glioblastoma multiform or GBM).
#'  
#'  This workflow introduces the following Bioconductor packages:
#'    
#'  * AnnotationHub
#'  * ChIPSeeker
#'  * ComplexHeatmap
#'  * pathview
#'  * ELMER
#'  * GAIA
#'  * MINET
#'  * RTCGAtoolbox
#'  * TCGAbiolinks
#'  * AnnotationHub
#'  
#'  Keywords: Epigenomics, Genomics, Cancer, non-coding, TCGA, ENCODE, Roadmap, Bioinformatics.
#' @docType package
#' @name TCGAWorkflow
#' @import TCGAbiolinks SummarizedExperiment ELMER AnnotationHub
#' @importFrom ChIPseeker tagMatrixList getTagMatrix
#' @importFrom gaia load_markers load_cnv runGAIA
#' @importFrom c3net c3net
#' @importFrom minet minet
#' @importFrom downloader download
#' @importFrom biomaRt getBM useMart 
#' @importFrom pathview pathview
#' @importFrom pbapply pbapply
#' @importFrom circlize circos.initializeWithIdeogram circos.genomicTrackPlotRegion circos.clear
#' @importFrom clusterProfiler bitr
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom ComplexHeatmap rowAnnotation HeatmapAnnotation Heatmap
#' @importFrom rGADEM GADEM consensus nOccurrences getPWM
#' @import  BSgenome.Hsapiens.UCSC.hg19 MotifDb   
#' @importFrom MotIV motifMatch
#' @importFrom motifStack plotMotifLogo
NULL
