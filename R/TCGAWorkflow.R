#' The aim of TCGAbiolinks is : 
#' i) facilitate the TCGA open-access data retrieval,
#' ii) prepare the data using the appropriate 
#' pre-processing strategies,
#' iii) provide the means to carry out different 
#' standard analyses and
#' iv) allow the user to download a specific version 
#' of the data and thus to easily reproduce earlier 
#' research results.
#' In more detail, the package provides multiple 
#' methods for analysis (e.g., differential expression 
#' analysis, identifying differentially methylated regions)
#' and methods for visualization (e.g., survival plots, 
#' volcano plots, starburst plots) in order to easily 
#' develop complete analysis pipelines.
#'
#' The functions you're likely to need from \pkg{TCGAbiolinks} is
#' \code{\link{GDCdownload}}, \code{\link{GDCquery}}.
#' Otherwise refer to the vignettes to see
#' how to format the documentation.
#'
#' @docType package
#' @name TCGAbiolinks
#' @import ELMER  downloader SummarizedExperiment  TCGAWorkflowExampleData
#' @import gaia ChIPseeker ComplexHeatmap
#' @import clusterProfiler minet c3net biomaRt pathview
#' @import BSgenome.Hsapiens.UCSC.hg19 motifStack
#' @import circlize pbapply GenomeInfoDb ggplot2 ggthemes parallel
#' @importFrom rGADEM GADEM nOccurrences getPWM consensus nMotifs
#' @importFrom MotIV viewAlignments motifMatch 
#' @importFrom AnnotationHub query AnnotationHub
#' @importFrom MotifDb MotifDb
NULL
