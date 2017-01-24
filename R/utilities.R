#' Creates a plot for GAIA ouptut (all significant aberrant regions.)
#' @description 
#' This function is a auxiliary function to visualize GAIA ouptut 
#' (all significant aberrant regions.)
#' @param calls A matrix with the following columns: Chromossome, Aberration Kind
#' Region Start, Region End, Region Size and score
#' @param threshold Score threshold (orange horizontal line in the plot)
#' @export
#' @examples 
#' call <- data.frame("Chromossome" = rep(9,100),
#'            "Aberration Kind" = rep(c(-2,-1,0,1,2),20),
#'            "Region Start [bp]" = 18259823:18259922,
#'            "Region End [bp]" = 18259823:18259922,
#'            "score"=rep(c(1,2,3,4),25))
#'  gaiaCNVplot(call,threshold = 0.01)  
#'  call <- data.frame("Chromossome" = rep(c(1,9),50),
#'            "Aberration Kind" = rep(c(-2,-1,0,1,2),20),
#'            "Region Start [bp]" = 18259823:18259922,
#'            "Region End [bp]" = 18259823:18259922,
#'            "score"=rep(c(1,2,3,4),25))
#'  gaiaCNVplot(call,threshold = 0.01)       
gaiaCNVplot <- function (calls,  threshold=0.01) {
  Calls <- calls[order(calls[,grep("start",colnames(calls),ignore.case = T)]),]
  Calls <- Calls[order(Calls[,grep("chr",colnames(calls),ignore.case = T)]),]
  rownames(Calls) <- NULL
  Chromo <- Calls[,grep("chr",colnames(calls),ignore.case = T)]
  Gains <- apply(Calls,1,function(x) ifelse(x[grep("aberration",colnames(calls),ignore.case = T)] == 1, x["score"], 0))
  Losses <- apply(Calls,1,function(x) ifelse(x[grep("aberration",colnames(calls),ignore.case = T)] == 0, x["score"], 0))
  plot(Gains, 
       ylim = c(-max(Calls[,"score"]+2), max(Calls[,"score"]+2)), 
       type = "h", 
       col = "red", 
       xlab = "Chromosome", 
       ylab = "Score", 
       xaxt = "n")
  points(-(Losses), type = "h", col = "blue")
  # Draw origin line
  abline(h = 0, cex = 4)
  # Draw threshold lines
  abline(h = -log10(threshold), col = "orange", cex = 4, main="test")
  abline(h = log10(threshold), col = "orange", cex = 4, main="test")
  
  uni.chr <- unique(Chromo)
  temp <- rep(0, length(uni.chr))
  for (i in 1:length(uni.chr)) {
    temp[i] <- max(which(uni.chr[i] == Chromo))
  }
  for (i in 1:length(temp)) {
    abline(v = temp[i], col = "black", lty = "dashed")
  }
  nChroms <- length(uni.chr)
  begin <- c()
  for (d in 1:nChroms) {
    chrom <- sum(Chromo == uni.chr[d])
    begin <- append(begin, chrom)
  }
  temp2 <- rep(0, nChroms)
  for (i in 1:nChroms) {
    if (i == 1) {
      temp2[1] <- (begin[1] * 0.5)
    }
    else if (i > 1) {
      temp2[i] <- temp[i - 1] + (begin[i] * 0.5)
    }
  }
  uni.chr[uni.chr==23] <- "X"
  uni.chr[uni.chr==24] <- "Y"
  for (i in 1:length(temp)) {
    axis(1, at = temp2[i], labels = uni.chr[i], cex.axis = 1)
  }
  legend(x=1,y=max(Calls[,"score"]+2), y.intersp=0.8, c("Amp"), pch=15, col=c("red"), text.font=3)
  legend(x=1,y=-max(Calls[,"score"]+0.5), y.intersp=0.8, c("Del"), pch=15, col=c("blue"), text.font=3)
}

#' Get a matrix of interactions of genes from biogrid
#' @description 
#' Using biogrid database, it will create a matrix of gene interations.
#' If columns A and row B has value 1, it means the gene A and gene B interatcs.
#' @param tmp.biogrid Biogrid table
#' @export
#' @param names.genes List of genes to filter from output. Default: consider all genes
#' @return A matrix with 1 for genes that interacts, 0 for no interaction.
get.adjacency.biogrid <- function(tmp.biogrid, names.genes = NULL){
  
  if(is.null(names.genes)){
    names.genes <- sort(union(unique(tmp.biogrid[,"Official.Symbol.Interactor.A"]), 
                              unique(tmp.biogrid[,"Official.Symbol.Interactor.B"])))
    ind <- seq(1,nrow(tmp.biogrid))
  }else{
    ind.A <- which(tmp.biogrid[,"Official.Symbol.Interactor.A"] %in% names.genes)
    ind.B <- which(tmp.biogrid[,"Official.Symbol.Interactor.B"] %in% names.genes)
    ind <- intersect(ind.A, ind.B)
  }
  
  mat.biogrid <- matrix(0, nrow=length(names.genes), ncol=length(names.genes), dimnames=list(names.genes, names.genes))
  
  for(i in ind){
    mat.biogrid[tmp.biogrid[i,"Official.Symbol.Interactor.A"], tmp.biogrid[i,"Official.Symbol.Interactor.B"]] <- mat.biogrid[tmp.biogrid[i,"Official.Symbol.Interactor.B"], tmp.biogrid[i,"Official.Symbol.Interactor.A"]] <- 1
  }
  diag(mat.biogrid) <- 0
  
  return(mat.biogrid)
}

#' Get GDC samples with both DNA methylation and Gene expression data
#' @description 
#' For a given TCGA project it gets the  samples (barcode) with both DNA methylation and Gene expression data
#' @param project A GDC project
#' @param n Number of samples to return. If NULL return all (default)
#' @usage matched_met_exp(project = "TCGA-ACC",n = 10)
#' @return A vector of barcodes
#' @importFrom TCGAbiolinks  GDCquery
#' @export
#' @examples 
#' # Get ACC samples with both DNA methylation and gene expression
#' samples <- matched_met_exp("TCGA-ACC")
matched_met_exp <- function(project, n = NULL){
  # get primary solid tumor samples: DNA methylation
  message("Download DNA methylation information")
  met450k <- GDCquery(project = project,
                      data.category = "DNA methylation",
                      platform = "Illumina Human Methylation 450",
                      legacy = TRUE, 
                      sample.type = c("Primary solid Tumor"))
  met450k.tp <-  met450k$results[[1]]$cases
  
  # get primary solid tumor samples: RNAseq
  message("Download gene expression information")
  exp <- GDCquery(project = project,
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = "Illumina HiSeq", 
                  file.type  = "results", 
                  sample.type = c("Primary solid Tumor"),
                  legacy = TRUE)
  exp.tp <-  exp$results[[1]]$cases
  print(exp.tp[1:10])
  # Get patients with samples in both platforms
  patients <- unique(substr(exp.tp,1,15)[substr(exp.tp,1,12) %in% substr(met450k.tp,1,12)] )
  if(!is.null(n)) patients <- patients[1:n] # get only n samples
  return(patients)
}