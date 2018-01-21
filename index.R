## ----setup, include=FALSE, echo=F, warning= F, message=F-------------------
knitr::opts_chunk$set(message = FALSE, 
                      warning = FALSE, 
                      error = FALSE, 
                      tidy = FALSE,
                      fig.align="center", 
                      dpi = 150, 
                      cache = FALSE,
                      progress=FALSE, 
                      quite = TRUE)

## ---- echo=FALSE, results="hide", warning=FALSE,message=FALSE--------------
suppressPackageStartupMessages({
  if("TCGAWorkflow" %in% installed.packages()[,1]) {
    library(TCGAWorkflow)
  } else {
    devtools::load_all(".")
  }
})

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  source("http://bioconductor.org/workflows.R")
#  workflowInstall("TCGAWorkflow")

## ---- eval=TRUE, message=FALSE, warning=FALSE, include=TRUE----------------
library(TCGAWorkflowData)
library(DT)

## ---- eval=TRUE, message=FALSE, warning=FALSE, include=FALSE---------------
library(TCGAbiolinks)

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  library(TCGAbiolinks)
#  # Obs: The data in the legacy database has been aligned to hg19
#  query.met.gbm <- GDCquery(project = "TCGA-GBM",
#                            legacy = TRUE,
#                            data.category = "DNA methylation",
#                            platform = "Illumina Human Methylation 450",
#                            barcode = c("TCGA-76-4926-01B-01D-1481-05", "TCGA-28-5211-01C-11D-1844-05"))
#  GDCdownload(query.met.gbm)
#  
#  met.gbm.450 <- GDCprepare(query = query.met.gbm,
#                            save = TRUE,
#                            save.filename = "gbmDNAmet450k.rda",
#                            summarizedExperiment = TRUE)
#  
#  query.met.lgg <- GDCquery(project = "TCGA-LGG",
#                            legacy = TRUE,
#                            data.category = "DNA methylation",
#                            platform = "Illumina Human Methylation 450",
#                            barcode = c("TCGA-HT-7879-01A-11D-2399-05", "TCGA-HT-8113-01A-11D-2399-05"))
#  GDCdownload(query.met.lgg)
#  met.lgg.450 <- GDCprepare(query = query.met.lgg,
#                            save = TRUE,
#                            save.filename = "lggDNAmet450k.rda",
#                            summarizedExperiment = TRUE)
#  met.gbm.lgg <- SummarizedExperiment::cbind(met.lgg.450, met.gbm.450)
#  
#  
#  query.exp.lgg <- GDCquery(project = "TCGA-LGG",
#                            legacy = TRUE,
#                            data.category = "Gene expression",
#                            data.type = "Gene expression quantification",
#                            platform = "Illumina HiSeq",
#                            file.type = "results",
#                            sample.type = "Primary solid Tumor")
#  GDCdownload(query.exp.lgg)
#  exp.lgg <- GDCprepare(query = query.exp.lgg, save = TRUE, save.filename = "lggExp.rda")
#  
#  query.exp.gbm <- GDCquery(project = "TCGA-GBM",
#                            legacy = TRUE,
#                            data.category = "Gene expression",
#                            data.type = "Gene expression quantification",
#                            platform = "Illumina HiSeq",
#                            file.type = "results",
#                            sample.type = "Primary solid Tumor")
#  GDCdownload(query.exp.gbm)
#  exp.gbm <- GDCprepare(query = query.exp.gbm, save = TRUE, save.filename = "gbmExp.rda")
#  exp.gbm.lgg <- SummarizedExperiment::cbind(exp.lgg, exp.gbm)

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  #-----------------------------------------------------------------------------
#  #                   Data.category: Copy number variation aligned to hg38
#  #-----------------------------------------------------------------------------
#  query <- GDCquery(project = "TCGA-ACC",
#                    data.category = "Copy Number Variation",
#                    data.type = "Copy Number Segment",
#                    barcode = c( "TCGA-OR-A5KU-01A-11D-A29H-01", "TCGA-OR-A5JK-01A-11D-A29H-01"))
#  GDCdownload(query)
#  data <- GDCprepare(query)
#  
#  query <- GDCquery("TCGA-ACC",
#                    "Copy Number Variation",
#                    data.type = "Masked Copy Number Segment",
#                    sample.type = c("Primary solid Tumor")) # see the barcodes with getResults(query)$cases
#  GDCdownload(query)
#  data <- GDCprepare(query)

## ---- eval=TRUE, include=TRUE----------------------------------------------
library(SummarizedExperiment)

# Load object from TCGAWorkflowData package
# THis object will be created in the further sections,
data(GBMIllumina_HiSeq) 

# get expression matrix
data <- assay(gbm.exp)
datatable(data[1:10,], 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = TRUE)

# get genes information
genes.info <- rowRanges(gbm.exp)
genes.info

# get sample information
sample.info <- colData(gbm.exp)
datatable(as.data.frame(sample.info), 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)


## ---- eval=TRUE, include=TRUE----------------------------------------------
# get indexed clinical patient data for GBM samples
gbm_clin <- GDCquery_clinic(project = "TCGA-GBM", type = "Clinical")

# get indexed clinical patient data for LGG samples
lgg_clin <- GDCquery_clinic(project = "TCGA-LGG", type = "Clinical")

# Bind the results, as the columns might not be the same,
# we will will plyr rbind.fill, to have all columns from both files
clinical <- plyr::rbind.fill(gbm_clin,lgg_clin)

## ----echo = TRUE, message = FALSE, warning = FALSE-------------------------
datatable(clinical[1:10,], options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-------------
# Fetch clinical data directly from the clinical XML files.
# if barcode is not set, it will consider all samples.
# We only set it to make the example faster
query <- GDCquery(project = "TCGA-GBM",
                  data.category = "Clinical",
                  barcode = c("TCGA-08-0516","TCGA-02-0317")) 
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")

## ----echo = TRUE, message = FALSE, warning = FALSE-------------------------
datatable(clinical, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-------------
clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug")

## ----echo = TRUE, message = FALSE, warning = FALSE-------------------------
datatable(clinical.drug, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-------------
clinical.radiation <- GDCprepare_clinic(query, clinical.info = "radiation")

## ----echo = TRUE, message = FALSE, warning = FALSE-------------------------
datatable(clinical.radiation, options = list(scrollX = TRUE,  keys = TRUE), rownames = FALSE)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-------------
clinical.admin <- GDCprepare_clinic(query, clinical.info = "admin")

## ----echo = TRUE, message = FALSE, warning = FALSE-------------------------
datatable(clinical.admin, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  LGGmut <- GDCquery_Maf(tumor = "LGG", pipelines = "mutect2")

## ----echo = TRUE, message = FALSE, warning = FALSE-------------------------
data(mafMutect2LGGGBM)
datatable(LGGmut[1:10,], options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

## ---- eval=TRUE, include=TRUE----------------------------------------------
gbm.subtypes <- TCGAquery_subtype(tumor = "gbm")

## ----echo = TRUE, message = FALSE, warning = FALSE-------------------------
datatable(gbm.subtypes[1:10,], options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  library(RTCGAToolbox)

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  # Get the last run dates
#  lastRunDate <- getFirehoseRunningDates()[1]
#  
#  # get DNA methylation data, RNAseq2 and clinical data for GBM
#  gbm.data <- getFirehoseData(dataset = "GBM",
#                              runDate = lastRunDate, gistic2_Date = getFirehoseAnalyzeDates(1),
#                              Methylation = FALSE, Clinic = TRUE,
#                              RNAseq2_Gene_Norm = FALSE, Mutation = TRUE,
#                              fileSizeLimit = 10000)
#  
#  gbm.mut <- getData(gbm.data,"Mutations")
#  gbm.clin <- getData(gbm.data,"Clinical")

## ---- eval=FALSE, message=FALSE, warning=FALSE, include=TRUE---------------
#  # Download GISTIC results
#  lastanalyzedate <- getFirehoseAnalyzeDates(1)
#  gistic <- getFirehoseData("GBM",gistic2_Date = lastanalyzedate)
#  
#  # get GISTIC results
#  gistic.allbygene <- getData(gistic,type = "GISTIC", CN = "All")
#  gistic.thresholedbygene <- getData(gistic,type = "GISTIC", CN = "Thresholed")

## ---- eval=TRUE, message=FALSE, warning=FALSE, include=TRUE----------------
data(GBMGistic)
datatable(gistic.allbygene,
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
datatable(gistic.thresholedbygene,
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

## ----results='hide', eval=FALSE, message=FALSE, warning=FALSE, include=TRUE----
#  library(TCGAbiolinks)
#  # Select common CN technology available for GBM and LGG
#  #############################
#  ## CNV data pre-processing ##
#  #############################
#  query.gbm.nocnv <- GDCquery(project = "TCGA-GBM",
#                              data.category = "Copy number variation",
#                              legacy = TRUE,
#                              file.type = "nocnv_hg19.seg",
#                              sample.type = c("Primary solid Tumor"))
#  # to reduce time we will select only 20 samples
#  query.gbm.nocnv$results[[1]] <- query.gbm.nocnv$results[[1]][1:20,]
#  
#  GDCdownload(query.gbm.nocnv, files.per.chunk = 100)
#  
#  gbm.nocnv <- GDCprepare(query.gbm.nocnv, save = TRUE, save.filename = "GBMnocnvhg19.rda")

## ----results='hide', eval=FALSE, message=FALSE, warning=FALSE, include=TRUE----
#  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-==--=--==---=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-==--=--==--
#  # Retrieve probes meta file from broadinstitute website for hg19
#  # For hg38 analysis please take a look on:
#  # https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files
#  # File: SNP6 GRCh38 Liftover Probeset File for Copy Number Variation Analysis
#  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-==--=--==---=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-==--=--==--
#  gdac.root <- "ftp://ftp.broadinstitute.org/pub/GISTIC2.0/hg19_support/"
#  file <- paste0(gdac.root, "genome.info.6.0_hg19.na31_minus_frequent_nan_probes_sorted_2.1.txt")
#  if(!file.exists(basename(file))) downloader::download(file, basename(file))
#  markersMatrix <-  readr::read_tsv(basename(file), col_names = FALSE, col_types = "ccn", progress = FALSE)
#  save(markersMatrix, file = "markersMatrix.rda", compress = "xz")

## ---- echo=TRUE, message=FALSE,warning=FALSE, include=TRUE-----------------
cancer <- "GBM"
message(paste0("Starting ", cancer))

# get objects created above
data(GBMnocnvhg19)
data(markersMatrix)

# Add label (0 for loss, 1 for gain)
cnvMatrix <- cbind(cnvMatrix,Label=NA)
cnvMatrix[cnvMatrix[,"Segment_Mean"] < -0.3,"Label"] <- 0
cnvMatrix[cnvMatrix[,"Segment_Mean"] > 0.3,"Label"] <- 1
cnvMatrix <- cnvMatrix[!is.na(cnvMatrix$Label),]

# Remove "Segment_Mean" and change col.names
cnvMatrix <- cnvMatrix[,-6]
colnames(cnvMatrix) <- c("Sample.Name", "Chromosome", "Start", "End", "Num.of.Markers", "Aberration")

# Substitute Chromosomes "X" and "Y" with "23" and "24"
cnvMatrix[cnvMatrix$Chromosome == "X","Chromosome"] <- 23
cnvMatrix[cnvMatrix$Chromosome == "Y","Chromosome"] <- 24
cnvMatrix$Chromosome <- as.integer(cnvMatrix$Chromosome)

# Recurrent CNV identification with GAIA
colnames(markersMatrix) <- c("Probe.Name", "Chromosome", "Start")
unique(markersMatrix$Chromosome)
markersMatrix[markersMatrix$Chromosome == "X","Chromosome"] <- 23
markersMatrix[markersMatrix$Chromosome == "Y","Chromosome"] <- 24
markersMatrix$Chromosome <- as.integer(markersMatrix$Chromosome)
markerID <- paste(markersMatrix$Chromosome,markersMatrix$Start, sep = ":")
# Removed duplicates
markersMatrix <- markersMatrix[!duplicated(markerID),]
# Filter markersMatrix for common CNV
markerID <- paste(markersMatrix$Chromosome,markersMatrix$Start, sep = ":")

file <- "ftp://ftp.broadinstitute.org/pub/GISTIC2.0/hg19_support/CNV.hg19.bypos.111213.txt"
if(!file.exists(basename(file))) downloader::download(file, basename(file))
commonCNV <- readr::read_tsv(basename(file), progress = FALSE)
commonID <- paste(commonCNV$Chromosome,commonCNV$Start, sep = ":")
markersMatrix_fil <- markersMatrix[!markerID %in% commonID,]

library(gaia)
set.seed(200)
markers_obj <- load_markers(as.data.frame(markersMatrix_fil))
nbsamples <- length(unique(cnvMatrix$Sample))
cnv_obj <- load_cnv(cnvMatrix, markers_obj, nbsamples)
suppressWarnings({
  results <- runGAIA(cnv_obj,
                     markers_obj,
                     output_file_name = paste0("GAIA_",cancer,"_flt.txt"),
                     aberrations = -1,  # -1 to all aberrations
                     chromosomes = 9, # -1 to all chromosomes
                     approximation = TRUE, # Set to TRUE to speed up the time requirements
                     num_iterations = 5000, # Reduced to speed up the time requirements
                     threshold = 0.25)
})
# Set q-value threshold
# Use a smalled value for your analysis. We set this as high values
# due to the small number of samples which did not reproduced
# results with smaller q-values
threshold <- 0.3

# Plot the results
RecCNV <- t(apply(results,1,as.numeric))
colnames(RecCNV) <- colnames(results)
RecCNV <- cbind(RecCNV, score = 0)
minval <- format(min(RecCNV[RecCNV[,"q-value"] != 0,"q-value"]), scientific = FALSE)
minval <- substring(minval,1, nchar(minval) - 1)
RecCNV[RecCNV[,"q-value"] == 0,"q-value"] <- as.numeric(minval)
RecCNV[,"score"] <- sapply(RecCNV[,"q-value"],function(x) -log10(as.numeric(x)))
RecCNV[RecCNV[,"q-value"] == as.numeric(minval),]

gaiaCNVplot(RecCNV,threshold)
save(results, RecCNV, threshold, file = paste0(cancer,"_CNV_results.rda"))

## ---- eval = TRUE, message=FALSE,warning=FALSE, include=FALSE--------------
library(GenomicRanges)

## ---- eval = FALSE, message=FALSE,warning=FALSE, include=TRUE--------------
#  library(GenomicRanges)
#  # Get gene information from GENCODE using biomart
#  genes <- TCGAbiolinks:::get.GRCh.bioMart(genome = "hg19")
#  genes <- genes[genes$external_gene_id != "" & genes$chromosome_name %in% c(1:22,"X","Y"),]
#  genes[genes$chromosome_name == "X", "chromosome_name"] <- 23
#  genes[genes$chromosome_name == "Y", "chromosome_name"] <- 24
#  genes$chromosome_name <- sapply(genes$chromosome_name,as.integer)
#  genes <- genes[order(genes$start_position),]
#  genes <- genes[order(genes$chromosome_name),]
#  genes <- genes[,c("external_gene_id", "chromosome_name", "start_position","end_position")]
#  colnames(genes) <- c("GeneSymbol","Chr","Start","End")
#  genes_GR <- makeGRangesFromDataFrame(genes,keep.extra.columns = TRUE)
#  save(genes_GR,genes,file = "genes_GR.rda", compress = "xz")

## ---- echo=TRUE, message=FALSE,warning=FALSE, include=TRUE-----------------
##############################
## Recurrent CNV annotation ## 
##############################
# Get gene information from GENCODE using biomart
data(genes_GR) # downloaded in the previous step (available in TCGAWorkflowData)

load(paste0(cancer,"_CNV_results.rda"))
sCNV <- RecCNV[RecCNV[,"q-value"] <= threshold,c(1:4,6)]
sCNV <- sCNV[order(sCNV[,3]),]
sCNV <- sCNV[order(sCNV[,1]),]
colnames(sCNV) <- c("Chr","Aberration","Start","End","q-value")
sCNV_GR <- makeGRangesFromDataFrame(sCNV,keep.extra.columns = TRUE)

hits <- findOverlaps(genes_GR, sCNV_GR, type = "within")
sCNV_ann <- cbind(sCNV[subjectHits(hits),],genes[queryHits(hits),])
AberrantRegion <- paste0(sCNV_ann[,1],":",sCNV_ann[,3],"-",sCNV_ann[,4])
GeneRegion <- paste0(sCNV_ann[,7],":",sCNV_ann[,8],"-",sCNV_ann[,9])
AmpDel_genes <- cbind(sCNV_ann[,c(6,2,5)],AberrantRegion,GeneRegion)
AmpDel_genes[AmpDel_genes[,2] == 0,2] <- "Del"
AmpDel_genes[AmpDel_genes[,2] == 1,2] <- "Amp"
rownames(AmpDel_genes) <- NULL

save(RecCNV, AmpDel_genes, file = paste0(cancer,"_CNV_results.rda"))

## ----results='asis', echo=FALSE, message=FALSE, warning=FALSE--------------
knitr::kable(head(AmpDel_genes), caption = "Chromosome 9 recurrent deleted genes in LGG")

## ----results='hide', echo=TRUE, message=FALSE, warning=FALSE, eval = FALSE----
#  library(maftools)
#  # Download Mutation Annotation Format (MAF) files
#  LGGmut <- GDCquery_Maf(tumor = "LGG", pipelines = "mutect2")
#  GBMmut <- GDCquery_Maf(tumor = "GBM", pipelines = "mutect2")
#  
#  # Merge them
#  mut <- plyr::rbind.fill(LGGmut, GBMmut)
#  save(maf,file ="mafMutect2LGGGBM.rda", compress = "xz")

## ----results='hide', echo=TRUE, message=FALSE, warning=FALSE---------------
library(maftools)
# recovering data from TCGAWorkflowData package.
data(mafMutect2LGGGBM)

# To prepare for maftools we will also include clinical data
# For a mutant vs WT survival analysis 
# get indexed clinical patient data for GBM samples
gbm_clin <- GDCquery_clinic(project = "TCGA-GBM", type = "Clinical")
# get indexed clinical patient data for LGG samples
lgg_clin <- GDCquery_clinic(project = "TCGA-LGG", type = "Clinical")
# Bind the results, as the columns might not be the same,
# we will will plyr rbind.fill, to have all columns from both files
clinical <- plyr::rbind.fill(gbm_clin,lgg_clin)
colnames(clinical)[1] <- "Tumor_Sample_Barcode"
clinical$Overall_Survival_Status <- 1
clinical$Overall_Survival_Status[which(clinical$vital_status != "dead")] <- 0
clinical$time <- clinical$days_to_death
clinical$time[is.na(clinical$days_to_death)] <- clinical$days_to_last_follow_up[is.na(clinical$days_to_death)]

# Create object to use in maftools
maf <- read.maf(maf = mut, clinicalData = clinical, isTCGA = T)

## ----echo=TRUE, message=FALSE, warning=FALSE,fig.width=10------------------
plotmafSummary(maf = maf,
               rmOutlier = TRUE,
               addStat = 'median',
               dashboard = TRUE)

## ----echo=TRUE, message=FALSE, warning=FALSE,fig.width=10------------------
gbm.subtypes <- TCGAquery_subtype(tumor = "gbm")
lgg.subtypes <- TCGAquery_subtype(tumor = "lgg")
subtype.info <- plyr::rbind.fill( lgg.subtypes, gbm.subtypes)
colnames(subtype.info)[1] <- "Tumor_Sample_Barcode" # required
oncoplot(maf = maf,
         top = 20,
         legendFontSize = 8,
         clinicalFeatures = c("Supervised.DNA.Methylation.Cluster","Original.Subtype"), # which columns should be plotted
         annotationDat = subtype.info) # data frame with metadata

## ----echo=TRUE, message=FALSE, warning=FALSE-------------------------------
plot <- mafSurvival(maf = maf,
                   genes = "IDH1",
                   fn = NULL,
                   time = 'time',
                   Status = 'Overall_Survival_Status',
                   isTCGA = TRUE)

## ----results='hide', eval = FALSE, echo=TRUE, message=FALSE,warning=FALSE----
#  LGGmut <- GDCquery_Maf(tumor = "LGG", pipelines = "mutect2")

## ----results='hide', echo=TRUE, message=FALSE,warning=FALSE----------------
###############################################
## Genomic aberration overview - Circos plot ## 
###############################################
# Retrieve curated mutations for selected cancer (e.g. "LGG") 
data(mafMutect2LGGGBM)
# Select only potentially damaging mutations
LGGmut <- LGGmut[LGGmut$Variant_Classification %in% c("Missense_Mutation",
                                                      "Nonsense_Mutation",
                                                      "Nonstop_Mutation",
                                                      "Frame_Shift_Del",
                                                      "Frame_Shift_Ins"),]
# Select recurrent mutations (identified in at least two samples)
mut.id <- paste0(LGGmut$Chromosome,":",LGGmut$Start_Position,"-",LGGmut$End_Position,"|",LGGmut$Reference_Allele)
mut <- cbind(mut.id, LGGmut)
# Prepare selected mutations data for circos plot
s.mut <- mut[mut$mut.id %in% unique(mut.id[duplicated(mut.id)]),]
s.mut <- s.mut[,c("Chromosome","Start_Position","End_Position","Variant_Classification","Hugo_Symbol")]
s.mut <- unique(s.mut)
typeNames <- unique(s.mut[,4])
type <- c(4:1)
names(type) <- typeNames[1:4]
Type <- type[s.mut[,4]]
s.mut <- cbind(s.mut,Type)
s.mut <- s.mut[,c(1:3,6,4,5)]

# Load recurrent CNV data for selected cancer (e.g. "LGG")
load("GBM_CNV_results.rda")
# Prepare selected sample CNV data for circos plot
s.cnv <- as.data.frame(RecCNV[RecCNV[,"q-value"] <= threshold,c(1:4,6)])
s.cnv <- s.cnv[,c(1,3,4,2)]
s.cnv[s.cnv$Chromosome == 23,"Chromosome"] <- "X"
s.cnv[s.cnv$Chromosome == 24,"Chromosome"] <- "Y"
Chromosome <- paste0("chr",s.cnv[,1])
s.cnv <- cbind(Chromosome, s.cnv[,-1])
s.cnv[,1] <- as.character(s.cnv[,1])
s.cnv[,4] <- as.character(s.cnv[,4])
s.cnv <- cbind(s.cnv,CNV=1)
colnames(s.cnv) <- c("Chromosome","Start_position","End_position","Aberration_Kind","CNV")

library(circlize)
# Draw genomic circos plot
par(mar=c(1,1,1,1), cex=1)
circos.initializeWithIdeogram()
# Add CNV results
colors <- c("forestgreen","firebrick")
names(colors)  <- c(0,1)
circos.genomicTrackPlotRegion(s.cnv,  ylim = c(0,1.2),
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, ytop.column = 2, ybottom = 0,
                                                   col = colors[value[[1]]], 
                                                   border="white")
                                cell.xlim = get.cell.meta.data("cell.xlim")
                                circos.lines(cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                              })
# Add mutation results
colors <- c("blue","green","red","gold")
names(colors)  <- typeNames[1:4]
circos.genomicTrackPlotRegion(s.mut, ylim = c(1.2,4.2),
                              panel.fun = function(region, value, ...) {
                                circos.genomicPoints(region, value, cex = 0.8, pch = 16, col = colors[value[[2]]], ...)
                              })

circos.clear()

legend(-0.2, 0.2, bty="n", y.intersp=1, c("Amp","Del"), pch=15, 
       col=c("firebrick","forestgreen"), title="CNVs", text.font=1, cex=0.4, title.adj=0)
legend(-0.2, 0, bty="n", y.intersp=1, names(colors), pch=16, 
       col=colors, title="Mutations", text.font=1, cex=0.4, title.adj=0)

## ----results='asis', echo=TRUE, message=FALSE,warning=FALSE----------------
par(mar=c(1,1,1,1),cex=1.5)
circos.par("start.degree" = 90, canvas.xlim = c(0, 1), canvas.ylim = c(0, 1), 
           gap.degree = 270, cell.padding = c(0, 0, 0, 0), track.margin = c(0.005, 0.005))
circos.initializeWithIdeogram(chromosome.index = "chr17")
circos.par(cell.padding = c(0, 0, 0, 0))
# Add CNV results
colors <- c("forestgreen","firebrick")
names(colors)  <- c(0,1)
circos.genomicTrackPlotRegion(s.cnv,  ylim = c(0,1.2),
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, ytop.column = 2, ybottom = 0,
                                                   col = colors[value[[1]]], 
                                                   border="white")
                                cell.xlim = get.cell.meta.data("cell.xlim")
                                circos.lines(cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                              })

# Add mutation results representing single genes
genes.mut <- paste0(s.mut$Hugo_Symbol,"-",s.mut$Type)
s.mutt <- cbind(s.mut,genes.mut)
n.mut <- table(genes.mut)
idx <- !duplicated(s.mutt$genes.mut)
s.mutt <- s.mutt[idx,]
s.mutt <- cbind(s.mutt,num=n.mut[s.mutt$genes.mut])
genes.num <- paste0(s.mutt$Hugo_Symbol," (",s.mutt$num.Freq,")")
s.mutt <- cbind(s.mutt[,-c(6:8)],genes.num)
s.mutt[,6] <- as.character(s.mutt[,6])
s.mutt[,4] <- s.mutt[,4]/2
s.mutt$num.Freq <- NULL
colors <- c("blue","green","red","gold")
names(colors)  <- typeNames[1:4]
circos.genomicTrackPlotRegion(s.mutt, ylim = c(0.3,2.2), track.height = 0.05,
                              panel.fun = function(region, value, ...) {
                                circos.genomicPoints(region, value, cex = 0.4, pch = 16, col = colors[value[[2]]], ...)
                              })

circos.genomicTrackPlotRegion(s.mutt, ylim = c(0, 1), track.height = 0.1, bg.border = NA)
i_track = get.cell.meta.data("track.index")

circos.genomicTrackPlotRegion(s.mutt, ylim = c(0,1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicText(region, value, 
                                                   y = 1, 
                                                   labels.column = 3,
                                                   col = colors[value[[2]]],
                                                   facing = "clockwise", adj = c(1, 0.5),
                                                   posTransform = posTransform.text, cex = 0.4, niceFacing = TRUE)
                              }, track.height = 0.1, bg.border = NA)

circos.genomicPosTransformLines(s.mutt,
                                posTransform = function(region, value)
                                  posTransform.text(region, 
                                                    y = 1, 
                                                    labels = value[[3]],
                                                    cex = 0.4, track.index = i_track+1),
                                direction = "inside", track.index = i_track)

circos.clear()

legend(0.25, 0.2, bty="n", y.intersp=1, c("Amp","Del"), pch=15, 
       col=c("firebrick","forestgreen"), title="CNVs", text.font=1, cex=0.4, title.adj=0)
legend(0, 0.2, bty="n", y.intersp=1, names(colors), pch=16, 
       col=colors, title="Mutations", text.font=1, cex=0.4, title.adj=0)

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  query <- GDCquery(project = "TCGA-GBM",
#                    data.category = "Gene expression",
#                    data.type = "Gene expression quantification",
#                    platform = "Illumina HiSeq",
#                    file.type  = "results",
#                    sample.type = c("Primary solid Tumor"),
#                    legacy = TRUE)
#  # We will use only 20 samples to make the example faster
#  query$results[[1]] <-  query$results[[1]][1:20,]
#  GDCdownload(query)
#  gbm.exp <- GDCprepare(query,
#                        save = TRUE,
#                        summarizedExperiment = TRUE,
#                        save.filename = "GBMIllumina_HiSeq.rda")
#  
#  query <- GDCquery(project = "TCGA-LGG",
#                    data.category = "Gene expression",
#                    data.type = "Gene expression quantification",
#                    platform = "Illumina HiSeq",
#                    file.type  = "results",
#                    sample.type = c("Primary solid Tumor"),
#                    legacy = TRUE)
#  # We will use only 20 samples to make the example faster
#  query$results[[1]] <-  query$results[[1]][1:20,]
#  GDCdownload(query)
#  lgg.exp <- GDCprepare(query,
#                        save = TRUE,
#                        summarizedExperiment = TRUE,
#                        save.filename = "LGGIllumina_HiSeq.rda")

## ----results='asis', echo=TRUE, message=FALSE,warning=FALSE----------------
data("LGGIllumina_HiSeq")
data("GBMIllumina_HiSeq")

dataPrep_LGG <- TCGAanalyze_Preprocessing(object = lgg.exp,
                                          cor.cut = 0.6,    
                                          datatype = "raw_count",
                                          filename = "LGG_IlluminaHiSeq_RNASeqV2.png")

dataPrep_GBM <- TCGAanalyze_Preprocessing(object = gbm.exp,
                                          cor.cut = 0.6, 
                                          datatype = "raw_count",
                                          filename = "GBM_IlluminaHiSeq_RNASeqV2.png")

dataNorm <- TCGAanalyze_Normalization(tabDF = cbind(dataPrep_LGG, dataPrep_GBM),
                                      geneInfo = TCGAbiolinks::geneInfo,
                                      method = "gcContent") #18323   672

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile",
                                  qnt.cut =  0.25)  # 13742   672

save(dataFilt, file = paste0("LGG_GBM_Norm_IlluminaHiSeq.rda"))

dataFiltLGG <- subset(dataFilt, select = substr(colnames(dataFilt),1,12) %in% lgg_clin$bcr_patient_barcode)
dataFiltGBM <- subset(dataFilt, select = substr(colnames(dataFilt),1,12) %in% gbm_clin$bcr_patient_barcode)

dataDEGs <- TCGAanalyze_DEA(mat1 = dataFiltLGG,
                            mat2 = dataFiltGBM,
                            Cond1type = "LGG",
                            Cond2type = "GBM",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")


## ----results='asis', echo=TRUE, message=FALSE,warning=FALSE----------------
# Number of differentially expressed genes (DEG)
nrow(dataDEGs)

## ----results='hide', echo=TRUE, message=FALSE,warning=FALSE,fig.height=10, fig.width=10----
#-------------------  4.2 EA: enrichment analysis             --------------------
ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes LGG Vs GBM", RegulonList = rownames(dataDEGs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
                        filename = NULL,
                        GOBPTab = ansEA$ResBP, 
                        nRGTab = rownames(dataDEGs),
                        nBar = 20)

## ----results='hide', echo=TRUE, message=FALSE,warning=FALSE,fig.height=10, fig.width=10----
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
                        filename = NULL,
                        GOCCTab = ansEA$ResCC,
                        nRGTab = rownames(dataDEGs),
                        nBar = 20)

## ----results='hide', echo=TRUE, message=FALSE,warning=FALSE,fig.height=10, fig.width=10----
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
                        filename = NULL,
                        GOMFTab = ansEA$ResMF, 
                        nRGTab = rownames(dataDEGs),
                        nBar = 20)

## ----results='hide', echo=TRUE, message=FALSE,warning=FALSE,fig.height=12, fig.width=15----
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
                        filename = NULL,
                        PathTab = ansEA$ResPat,
                        nRGTab = rownames(dataDEGs),
                        nBar = 20)

## ----results='asis', echo=TRUE, message=FALSE,warning=FALSE----------------
library(SummarizedExperiment)
GenelistComplete <- rownames(assay(lgg.exp,1))

# DEGs TopTable
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"LGG","GBM",
                                          dataFilt[,colnames(dataFiltLGG)],
                                          dataFilt[,colnames(dataFiltGBM)])

dataDEGsFiltLevel$GeneID <- 0

library(clusterProfiler)
# Converting Gene symbol to geneID
eg = as.data.frame(bitr(dataDEGsFiltLevel$mRNA,
                        fromType="SYMBOL",
                        toType="ENTREZID",
                        OrgDb="org.Hs.eg.db"))
eg <- eg[!duplicated(eg$SYMBOL),]

dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]

dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]

# table(eg$SYMBOL == dataDEGsFiltLevel$mRNA) should be TRUE
all(eg$SYMBOL == dataDEGsFiltLevel$mRNA)
dataDEGsFiltLevel$GeneID <- eg$ENTREZID

dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
library(pathview)
# pathway.id: hsa05214 is the glioma pathway
# limit: sets the limit for gene expression legend and color
hsa05214 <- pathview::pathview(gene.data  = genelistDEGs,
                               pathway.id = "hsa05214",
                               species    = "hsa",
                               limit = list(gene=as.integer(max(abs(genelistDEGs)))))


## ---- eval=FALSE, include=TRUE---------------------------------------------
#  ### read biogrid info (available in TCGAWorkflowData as "biogrid")
#  ### Check last version in https://thebiogrid.org/download.php
#  file <- "http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.4.146/BIOGRID-ALL-3.4.146.tab2.zip"
#  if(!file.exists(gsub("zip","txt",basename(file)))){
#    downloader::download(file,basename(file))
#    unzip(basename(file),junkpaths =TRUE)
#  }
#  
#  tmp.biogrid <- read.csv(gsub("zip","txt",basename(file)), header=TRUE, sep="\t", stringsAsFactors=FALSE)

## ---- eval=TRUE, include=TRUE----------------------------------------------
### plot details (colors & symbols)
mycols <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628')

### load network inference libraries
library(minet)
library(c3net)

### deferentially identified genes using TCGAbiolinks
# we will use only a subset (first 50 genes) of it to make the example faster
names.genes.de <- rownames(dataDEGs)[1:30]

data(biogrid)
net.biogrid.de <- getAdjacencyBiogrid(tmp.biogrid, names.genes.de)

mydata <- dataFiltLGG[names.genes.de, ]

### infer networks
t.mydata <- t(mydata)
net.aracne <- minet(t.mydata, method = "aracne")
net.mrnet <- minet(t.mydata)
net.clr <- minet(t.mydata, method = "clr")
net.c3net <- c3net(mydata)

### validate compared to biogrid network
tmp.val <- list(validate(net.aracne, net.biogrid.de), 
                validate(net.mrnet, net.biogrid.de),
                validate(net.clr, net.biogrid.de), 
                validate(net.c3net, net.biogrid.de))

### plot roc and compute auc for the different networks
dev1 <- show.roc(tmp.val[[1]],cex=0.3,col=mycols[1],type="l")
res.auc <- auc.roc(tmp.val[[1]])
for(count in 2:length(tmp.val)){
  show.roc(tmp.val[[count]],device=dev1,cex=0.3,col=mycols[count],type="l")
  res.auc <- c(res.auc, auc.roc(tmp.val[[count]]))
}

legend("bottomright", legend=paste(c("aracne","mrnet","clr","c3net"), signif(res.auc,4), sep=": "),
       col=mycols[1:length(tmp.val)],lty=1, bty="n" )
# Please, uncomment this line to produce the pdf files.
# dev.copy2pdf(width=8,height=8,device = dev1, file = paste0("roc_biogrid_",cancertype,".pdf"))

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  #----------------------------
#  # Obtaining DNA methylation
#  #----------------------------
#  # Samples
#  lgg.samples <- matchedMetExp("TCGA-LGG", n = 10)
#  gbm.samples <- matchedMetExp("TCGA-GBM", n = 10)
#  samples <- c(lgg.samples,gbm.samples)
#  
#  #-----------------------------------
#  # 1 - Methylation
#  # ----------------------------------
#  # For methylation it is quicker in this case to download the tar.gz file
#  # and get the samples we want instead of downloading files by files
#  query <- GDCquery(project = c("TCGA-LGG","TCGA-GBM"),
#                    data.category = "DNA methylation",
#                    platform = "Illumina Human Methylation 450",
#                    legacy = TRUE,
#                    barcode = samples)
#  GDCdownload(query)
#  met <- GDCprepare(query, save = FALSE)
#  
#  # We will use only chr9 to make the example faster
#  met <- subset(met,subset = as.character(seqnames(met)) %in% c("chr9"))
#  # This data is avaliable in the package (object elmerExample)

## ----echo=TRUE, message=FALSE,warning=FALSE, fig.height=8,fig.width=8------
data(elmerExample)
#----------------------------
# Mean methylation
#----------------------------
# Plot a barplot for the groups in the disease column in the
# summarizedExperiment object

# remove probes with NA (similar to na.omit)
met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0))

TCGAvisualize_meanMethylation(met,
                              groupCol = "project_id",
                              group.legend  = "Groups",
                              filename = NULL,
                              print.pvalue = TRUE)

## ----results='hide', echo=TRUE, message=FALSE,warning=FALSE----------------
#------- Searching for differentially methylated CpG sites     ----------
met <- TCGAanalyze_DMR(met,
                       groupCol = "project_id", # a column in the colData matrix
                       group1 = "TCGA-GBM", # a type of the disease type column
                       group2 = "TCGA-LGG", # a type of the disease column
                       p.cut = 0.05,
                       diffmean.cut = 0.15,
                       save = FALSE,
                       legend = "State",
                       plot.filename = "LGG_GBM_metvolcano.png",
                       cores = 1 # if set to 1 there will be a progress bar
)

## ----results='hide', echo=TRUE, message=FALSE,warning=FALSE----------------
#--------------------------
# DNA Methylation heatmap
#-------------------------
library(ComplexHeatmap)
clinical <- plyr::rbind.fill(gbm_clin,lgg_clin)

# get the probes that are Hypermethylated or Hypomethylated
# met is the same object of the section 'DNA methylation analysis'
status.col <- "status.TCGA.GBM.TCGA.LGG"
sig.met <- met[values(met)[,status.col] %in% c("Hypermethylated","Hypomethylated"),]


# top annotation, which sampples are LGG and GBM
# We will add clinical data as annotation of the samples
# we will sort the clinical data to have the same order of the DNA methylation matrix
clinical.order <- clinical[match(substr(colnames(sig.met),1,12),clinical$bcr_patient_barcode),]
ta = HeatmapAnnotation(df = clinical.order[,c("disease","gender","vital_status","race")],
                       col = list(disease = c("LGG" = "grey", "GBM" = "black"),
                                  gender = c("male" = "blue","female" = "pink")))

# row annotation: add the status for LGG in relation to GBM
# For exmaple: status.gbm.lgg Hypomethyated means that the
# mean DNA methylation of probes for lgg are hypomethylated
# compared to GBM ones.
ra = rowAnnotation(df = values(sig.met)[status.col],
                   col = list("status.TCGA.GBM.TCGA.LGG" = 
                                c("Hypomethylated" = "orange",
                                  "Hypermethylated" = "darkgreen")),
                   width = unit(1, "cm"))

heatmap  <- Heatmap(assay(sig.met),
                    name = "DNA methylation",
                    col = matlab::jet.colors(200),
                    show_row_names = FALSE,
                    cluster_rows = TRUE,
                    cluster_columns = FALSE,
                    show_column_names = FALSE,
                    bottom_annotation = ta,
                    column_title = "DNA Methylation") 
# Save to pdf
png("heatmap.png",width = 600, height = 400)
draw(heatmap, annotation_legend_side =  "bottom")
dev.off()

## ----results='asis', echo=TRUE, message=FALSE,warning=FALSE----------------
library(rGADEM)
library(BSgenome.Hsapiens.UCSC.hg19)
library(motifStack)

probes <- rowRanges(met)[values(met)[,status.col] %in% c("Hypermethylated" ,"Hypomethylated") ,]
# Get hypo/hyper methylated probes and make a 200bp window 
# surrounding each probe.
sequence <- RangedData(space = as.character(seqnames(probes)),
                       IRanges(start = start(ranges(probes)) - 100,
                               end = start(ranges(probes)) + 100), strand = "*")
#look for motifs
gadem <- GADEM(sequence, verbose = FALSE, genome = Hsapiens)

# How many motifs were found?
nMotifs(gadem)

# get the number of occurrences
nOccurrences(gadem)

# view all sequences consensus
consensus(gadem)

# Print motif
pwm <- getPWM(gadem)
pfm  <- new("pfm",mat=pwm[[1]],name="Novel Site 1")
plotMotifLogo(pfm)

# Number of instances of motif 1?
length(gadem@motifList[[1]]@alignList)


## ----results='asis', echo=TRUE,message=FALSE,warning=FALSE-----------------
library(MotIV)
# motif matching analysis
analysis.jaspar <- motifMatch(pwm)
summary(analysis.jaspar)
plot(analysis.jaspar, ncol=2, top=5, rev=FALSE, main="", bysim=TRUE, cex=0.4)

#  visualize the quality of the results looking into the alignments
# E-value gives an estimation of the match.
alignment <- viewAlignments(analysis.jaspar)
print(alignment[[1]])

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  #------------------- Starburst plot ------------------------------
#  starburst <- TCGAvisualize_starburst(met,    # DNA methylation with results
#                                       dataDEGs,    # DEG results
#                                       group1 = "Glioblastoma Multiforme",
#                                       group2 = "Brain Lower Grade Glioma",
#                                       filename = "starburst.png",
#                                       genome = "hg19",
#                                       met.platform = "450K",
#                                       met.p.cut = 0.05,
#                                       exp.p.cut = 10^-2,
#                                       diffmean.cut = 0.15,
#                                       logFC.cut = 1,
#                                       width = 15,height = 10,
#                                       names = TRUE)

## ----results='asis', echo=FALSE, message=FALSE, warning=FALSE--------------
#------------------- Starburst plot ------------------------------
starburst <- TCGAvisualize_starburst(met,    # DNA methylation with results
                                     dataDEGs,    # DEG results
                                     group1 = "TCGA-GBM", 
                                     group2 = "TCGA-LGG", 
                                     filename = "starburst.png",
                                     genome = "hg19",
                                     met.platform = "450K",
                                     met.p.cut = 0.05,
                                     exp.p.cut = 10^-2,
                                     diffmean.cut = 0.15,
                                     logFC.cut = 1,
                                     width = 15,height = 10,
                                     names = TRUE)

## ----table4, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'-----
tabl <- "
|                  Histone marks                 |                                                   Role                                                  |
|:----------------------------------------------:|:-------------------------------------------------------------------------------------------------------:|
| Histone H3 lysine 4 trimethylation (H3K4me3)   | Promoter regions [@heintzman2007distinct,@bernstein2005genomic]                                      |
| Histone H3 lysine 4 monomethylation (H3K4me1)  | Enhancer regions [@heintzman2007distinct]                                                           |
| Histone H3 lysine 36 trimethylation (H3K36me3) | Transcribed regions                                                                                     |
| Histone H3 lysine 27 trimethylation (H3K27me3) | Polycomb repression [@bonasio2010molecular]                                                         |
| Histone H3 lysine 9 trimethylation (H3K9me3)   | Heterochromatin regions  [@peters2003partitioning]                                                  |
| Histone H3 acetylated at lysine 27 (H3K27ac)   | Increase activation of genomic elements [@heintzman2009histone,@rada2011unique,@creyghton2010histone] |
| Histone H3 lysine 9 acetylation  (H3K9ac)      | Transcriptional activation [@nishida2006histone]                                                    |
"
cat(tabl) 

## ----table5, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'-----
tabl <- "  
|                File                |                               Description                              |
|:----------------------------------:|:----------------------------------------------------------------------:|
| fc.signal.bigwig                   | Bigwig File containing  fold enrichment signal tracks                  |
| pval.signal.bigwig                 | Bigwig File containing -log10(p-value) signal tracks                   |
| hotspot.fdr0.01.broad.bed.gz       | Broad domains on enrichment for  DNase-seq for consolidated epigenomes |
| hotspot.broad.bed.gz               | Broad domains on enrichment for DNase-seq for consolidated epigenomes  |
| broadPeak.gz                       | Broad ChIP-seq peaks for consolidated  epigenomes                      |
| gappedPeak.gz                      | Gapped ChIP-seq peaks for consolidated   epigenomes                    |
| narrowPeak.gz                      | Narrow ChIP-seq peaks for consolidated epigenomes                      |
| hotspot.fdr0.01.peaks.bed.gz       | Narrow DNasePeaks for   consolidated epigenomes                        |
| hotspot.all.peaks.bed.gz           | Narrow DNasePeaks for  consolidated epigenomes                         |
| .macs2.narrowPeak.gz               | Narrow DNasePeaks for consolidated epigenomes                          |
| coreMarks_mnemonics.bed.gz        | 15 state chromatin segmentations                                       |
| mCRF_FractionalMethylation.bigwig | MeDIP/MRE(mCRF) fractional methylation calls                           |
| RRBS_FractionalMethylation.bigwig | RRBS fractional methylation calls                                      |
| WGBS_FractionalMethylation.bigwig | Whole genome bisulphite fractional methylation calls                   |
"
cat(tabl) 

## ----results='hide', eval=TRUE, echo=TRUE, message=FALSE,warning=FALSE-----
library(ChIPseeker)
library(pbapply)
library(ggplot2)

## ----results='hide', eval=FALSE, echo=TRUE, message=FALSE,warning=FALSE----
#  #------------------ Working with ChipSeq data ---------------
#  # Step 1: download histone marks for a brain and non-brain samples.
#  #------------------------------------------------------------
#  # loading annotation hub database
#  library(AnnotationHub)
#  ah = AnnotationHub()
#  
#  # Searching for brain consolidated epigenomes in the roadmap database
#  bpChipEpi_brain <- query(ah , c("EpigenomeRoadMap", "narrowPeak", "chip", "consolidated","brain","E068"))
#  # Get chip-seq data
#  histone.marks <- pblapply(names(bpChipEpi_brain), function(x) {ah[[x]]})
#  names(histone.marks) <- names(bpChipEpi_brain)
#  # OBS: histone.marks is available in TCGAWorkflowData package

## ----results='hide', echo=TRUE, message=FALSE,warning=FALSE----------------
data(histoneMarks)
# Create a GR object based on the hypo/hypermethylated probes.
probes <- keepStandardChromosomes(rowRanges(met)[values(met)[,status.col] %in% c("Hypermethylated","Hypomethylated"),])
# Defining a window of 3kbp - 3kbp_probe_3kbp
ranges(probes) <- IRanges(start = c(start(ranges(probes)) - 3000), end = c(start(ranges(probes)) + 3000))

### Profile of ChIP peaks binding to TSS regions
# First of all, to calculate the profile of ChIP peaks binding to TSS regions, we should
# prepare the TSS regions, which are defined as the flanking sequence of the TSS sites.
# Then align the peaks that are mapping to these regions and generate the tagMatrix.
tagMatrixList <- pbapply::pblapply(histone.marks, function(x) {
  getTagMatrix(keepStandardChromosomes(x), windows = probes, weightCol = "score")
})
# change names retrieved with the following command: basename(bpChipEpi_brain$title)
names(tagMatrixList) <- c("H3K4me1","H3K4me3", "H3K9ac", "H3K9me3", "H3K27ac",  "H3K27me3", "H3K36me3")

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  tagHeatmap(tagMatrixList, xlim=c(-3000, 3000),color = NULL)

## ----results='asis', echo=FALSE, fig.width = 7, message=FALSE,warning=FALSE----
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000),color = NULL)

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  p <- plotAvgProf(tagMatrixList, xlim = c(-3000,3000), xlab = "Genomic Region (5'->3', centered on CpG)")
#  # We are centreing in the CpG instead of the TSS. So we'll change the labels manually
#  p <- p + scale_x_continuous(breaks=c(-3000,-1500,0,1500,3000),labels=c(-3000,-1500,"CpG",1500,3000))
#  library(ggthemes)
#  p + theme_few() + scale_colour_few(name="Histone marks") +  guides(colour = guide_legend(override.aes = list(size=4)))

## ----results='asis', echo=FALSE, message=FALSE,warning=FALSE---------------
p <- plotAvgProf(tagMatrixList, xlim = c(-3000,3000), xlab = "Genomic Region (5'->3', centered on CpG)")
# We are centreing in the CpG instead of the TSS. So we'll change the labels manually
if (requireNamespace("ggplot2", quietly = TRUE)) 
  p <- p + ggplot2::scale_x_continuous(breaks=c(-3000,-1500,0,1500,3000),
                                       labels=c(-3000,-1500,"CpG",1500,3000))

if (requireNamespace("ggthemes", quietly = TRUE)) 
  p + ggthemes::theme_few() + 
  ggthemes::scale_colour_few(name="Histone marks") + 
  ggplot2::guides(colour = guide_legend(override.aes = list(size=4)))

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  #----------- 8.3 Identification of Regulatory Enhancers   -------
#  library(TCGAbiolinks)
#  # Samples: primary solid tumor w/ DNA methylation and gene expression
#  lgg.samples <- matchedMetExp("TCGA-LGG", n = 10)
#  gbm.samples <- matchedMetExp("TCGA-GBM", n = 10)
#  samples <- c(lgg.samples,gbm.samples)
#  
#  #-----------------------------------
#  # 1 - Methylation
#  # ----------------------------------
#  query.met <- GDCquery(project = c("TCGA-LGG","TCGA-GBM"),
#                        data.category = "DNA methylation",
#                        platform = "Illumina Human Methylation 450",
#                        legacy = TRUE,
#                        barcode = samples)
#  GDCdownload(query.met)
#  met <- GDCprepare(query.met, save = FALSE)
#  met <- subset(met,subset = as.character(GenomicRanges::seqnames(met)) %in% c("chr9"))
#  
#  #-----------------------------------
#  # 2 - Expression
#  # ----------------------------------
#  query.exp <- GDCquery(project = c("TCGA-LGG","TCGA-GBM"),
#                        data.category = "Gene expression",
#                        data.type = "Gene expression quantification",
#                        platform = "Illumina HiSeq",
#                        file.type  = "results",
#                        legacy = TRUE,
#                        barcode =  samples)
#  GDCdownload(query.exp)
#  exp <- GDCprepare(query.exp, save = FALSE)
#  save(exp, met, gbm.samples, lgg.samples, file = "elmer.example.rda", compress = "xz")

## ---- message=FALSE,warning=FALSE------------------------------------------
library(ELMER)
library(MultiAssayExperiment)
library(GenomicRanges)
distal.probes <- get.feature.probe(genome = "hg19", 
                                   met.platform = "450K")

# Recover the data created in the last step
data(elmerExample)
rownames(exp) <- values(exp)$ensembl_gene_id
mae <- createMAE(exp = exp, 
                 met = met,
                 save = TRUE,
                 linearize.exp = TRUE,
                 save.filename = "mae.rda",
                 filter.probes = distal.probes,
                 met.platform = "450K",
                 genome = "hg19",
                 TCGA = TRUE)
mae

## ---- warning=FALSE, results="hide"----------------------------------------
cores <- 1 # you can use more cores if you want
group.col <- "project_id"
group1 <- "TCGA-GBM"
group2 <- "TCGA-LGG"

# Available directions are hypo and hyper, we will use only hypo
# due to speed constraint
direction <- "hypo"
dir.out <- paste0("elmer/",direction)
dir.create(dir.out, showWarnings = FALSE, recursive = TRUE)

## ---- warning=FALSE, results="hide"----------------------------------------
#--------------------------------------
# STEP 3: Analysis                     |
#--------------------------------------
# Step 3.1: Get diff methylated probes |
#--------------------------------------
message("Get diff methylated probes")
Sig.probes <- get.diff.meth(data = mae, 
                            group.col = group.col,
                            group1 = group1,
                            group2 =  group2,
                            minSubgroupFrac = 1.0, # Use all samples
                            sig.dif = 0.2, # defualt is 0.3
                            diff.dir = direction, # Search for hypomethylated probes in group 1
                            cores = cores, 
                            dir.out = dir.out, 
                            pvalue = 0.1)
getMethdiff <- readr::read_csv(paste0(dir.out,"/getMethdiff.",direction,".probes.significant.csv")  )

## ---- warning=FALSE--------------------------------------------------------
datatable(getMethdiff[1:10,], 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = TRUE)

## ---- warning=FALSE, results="hide"----------------------------------------
#-------------------------------------------------------------
# Step 3.2: Identify significant probe-gene pairs            |
#-------------------------------------------------------------
# Collect nearby 20 genes for Sig.probes
message("Get nearby genes")
nearGenes <- GetNearGenes(data = mae,
                          numFlankingGenes = 4, # default is 20 genes
                          probes = Sig.probes$probe)

message("Get anti correlated probes-genes")
pair <- get.pair(data = mae,
                 group.col = group.col,
                 group1 = group1,
                 group2 =  group2,
                 nearGenes = nearGenes,
                 mode = "supervised",
                 minSubgroupFrac = 1, # % of samples to use in to create groups U/M
                 permu.dir = paste0(dir.out,"/permu"),
                 permu.size = 5, # Please set to 100000 to get significant results
                 raw.pvalue = 0.1,   # defualt is 0.001
                 Pe = 0.5, # Please set to 0.001 to get significant results
                 filter.probes = TRUE, # See preAssociationProbeFiltering function
                 filter.percentage = 0.05,
                 filter.portion = 0.3,
                 dir.out = dir.out,
                 diff.dir = direction,
                 cores = cores,
                 label = direction)

## ---- warning=FALSE--------------------------------------------------------
nearGenes[[1]]
datatable(pair[1:10,], 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = TRUE)

## ---- warning=FALSE, results="hide"----------------------------------------
#-------------------------------------------------------------
# Step 3.3: Motif enrichment analysis on the selected probes |
#-------------------------------------------------------------
enriched.motif <- get.enriched.motif(data = mae,
                                     probes = pair$Probe, 
                                     dir.out = dir.out, 
                                     pvalue = 0.7, # default is 0.05
                                     label = direction,
                                     min.incidence = 10,
                                     lower.OR = 1.1)
# One of the output from the  previous function is a file with the motif, OR and Number of probes
# It will be used for plotting purposes
motif.enrichment <- read.csv(paste0(dir.out,"/getMotif.",direction, ".motif.enrichment.csv"))

## ---- warning=FALSE--------------------------------------------------------
head(enriched.motif[names(enriched.motif)[1]]) ## probes in the given set that have the first motif.
datatable(motif.enrichment, 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          filter = 'top',
          rownames = FALSE)

## ---- warning=FALSE, results="hide"----------------------------------------
#-------------------------------------------------------------
# Step 3.4: Identifying regulatory TFs                        |
#-------------------------------------------------------------
TF <- get.TFs(data = mae, 
              group.col = group.col,
              group1 = group1,
              group2 =  group2,
              mode = "supervised",
              enriched.motif = enriched.motif,
              dir.out = dir.out, 
              cores = cores, 
              diff.dir = direction,
              label = direction)
    
# One of the output from the previous function is a file with the raking of TF,
# for each motif. It will be used for plotting purposes
TF.meth.cor <- get(load(paste0(dir.out,"/getTF.",direction,".TFs.with.motif.pvalue.rda")))

## ---- warning=FALSE--------------------------------------------------------
datatable(TF, 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
datatable(TF.meth.cor[1:10,1:10], 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = TRUE)

## ---- echo=TRUE, message=FALSE, warnings=FALSE, results='asis',fig.height=8, fig.width=10----
heatmapPairs(data = mae, 
             group.col = group.col,
             group1 = group1, 
             group2 = group2, 
             annotation.col = c("gender"),
             pairs = pair,
             filename =  NULL)

## ---- echo=TRUE, message=FALSE, warnings=FALSE, results='asis',fig.height=8,fig.width=15----
motif.enrichment.plot(motif.enrichment = motif.enrichment, 
                      save = FALSE,
                      significant = list(lowerOR = 1.1)) # Filter motifs in the plot lowerOR > 1.3

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  grid:TF.rank.plot(motif.pvalue=TF.meth.cor, motif=TF$motif[1], save=FALSE)

## ----results='asis', echo=FALSE, message=FALSE,warning=FALSE,  fig.align='center',fig.height=8,fig.width=8----
gg <- TF.rank.plot(motif.pvalue=TF.meth.cor, motif=TF$motif[1], save=FALSE)
grid::grid.draw(gg[[1]])

## ----results='asis', echo=TRUE, message=FALSE, warning=FALSE---------------
png("TF.png",width = 800, height = 400)
scatter.plot(mae, 
             category = group.col, 
             save = FALSE, 
             lm_line = TRUE,
             byTF=list(TF=unlist(stringr::str_split(TF[1,"top_5percent_TFs"],";"))[1:4], 
                       probe=enriched.motif[[TF$motif[1]]]))
dev.off()

## --------------------------------------------------------------------------
pander::pander(sessionInfo(), compact = FALSE)

