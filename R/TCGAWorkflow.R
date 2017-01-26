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
#' @import ELMER AnnotationHub downloader SummarizedExperiment 
#' @import gaia ChIPseeker ComplexHeatmap
#' @import clusterProfiler minet c3net biomaRt pathview
#' @import BSgenome.Hsapiens.UCSC.hg19 MotifDb MotIV rGADEM motifStack
#' @import circlize pbapply GenomeInfoDb ggplot2 ggthemes parallel
NULL

#' A SummarizedExperiment containing
#' TCGA data: DNA methylation platform 450K chromossome 9 
#' for 10 LGG samples and 10 GBM samples
#' @docType data
#' @keywords internal
#' @name met
#' @format A SumarrizedExperiment with 9861 rows and 20 samples
NULL

#' A DNA methylation matrix for 10 GBM and 10 LGG samples prepared for
#' the creation of an ELMER object.
#' @docType data
#' @keywords internal
#' @name met.elmer
#' @format A matrix with 8066 rows and 20 samples
NULL

#' A gene expression matrix for 10 GBM and 10 LGG samples prepared for
#' the creation of an ELMER object.
#' @docType data
#' @keywords internal
#' @name exp.elmer
#' @format A matrix with 21022 rows and 20 samples
NULL

#' Identifiers for the 10 LGG samples in the ELMER objects
#' @docType data
#' @keywords internal
#' @name lgg.samples
#' @format A vector of 10 barcodes
NULL

#' Identifiers for the 10 GBM samples in the ELMER objects
#' @docType data
#' @keywords internal
#' @name gbm.samples
#' @format A vector of 10 barcodes
NULL

#' A gene expression matrix wih 20 LGG samples
#' @docType data
#' @keywords internal
#' @name lgg.exp
#' @format Gene expression: A SummariedExperiment object 
#' with 21022 rows and 20 columns
NULL

#' A gene expression matrix wih 20 GBM samples
#' @docType data
#' @keywords internal
#' @name gbm.exp
#' @format Gene expression: A SummariedExperiment object
#'  with 21022 rows and 20 columns
NULL

#' CNV data for 20 TCGA-GBM samples
#' @docType data
#' @keywords internal
#' @name cnvMatrix
#' @format A matrix with 3252 rows and 6 columns
NULL