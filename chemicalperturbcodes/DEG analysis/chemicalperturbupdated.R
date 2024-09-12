#readRDS file
endoc<-readRDS("/athena/chenlab/scratch/jjv4001/endocchemperturbnew.RDS")
library(Seurat)
library(DESeq2)
library(tidyverse)
#normalization, scaling and clustering
degAnalysis <- function(hormone1, hormone2) {
  #Choose hormone 1 and subset cells from endoc object to a separate object, create a column in metadata called cell_label to label cell 1 to cell n corresponding to the number of cells in that object
  H1<-subset(endoc, hormone== hormone1)
  num_cells_H1 <- ncol(H1)
  H1$cell_label <- paste0("cell", 1:num_cells_H1)
  #Choose hormone 2 and subset cells from endoc object to a separate object, create a column in metadata called cell_label to label cell 1 to cell n corresponding to the number of cells in that object
  H2<-subset(endoc, hormone== hormone2)
  num_cells_H2 <- ncol(H2)
  H2$cell_label <- paste0("cell", 1:num_cells_H2)
  #merge objects to create one count matrix
  combined<-merge(H1,H2)
  #run Aggregate expression to perform pseudobulk DE seq, this aggregates counts for all genes according to the groups supplied
  counts<-AggregateExpression(combined, group.by= c("cell_label","hormone"),slot= "counts", return.seurat = FALSE)
  counts<-counts$RNA
  #order columns in counts matrix according to names 
  suffix <- gsub(".*_", "", colnames(counts))
  sorted_colnames <- colnames(counts)[order(suffix)]
  counts <- counts[, sorted_colnames]
  #create colData corresponding to each cell and hormone to be distinguished (col1=cell1_hormone1, col2=hormone1 or col1=cell100_hormone2, col2=hormone2)
  colData<-data.frame(samples=colnames(counts))
  rownames(colData) <- NULL
  colData$hormone<-gsub('.*_','',colData$samples)
  colData <- colData %>% 
    column_to_rownames(var = 'samples')
  #runDESeq2 analysis
  dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = colData, design = ~hormone)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds<-DESeq(dds)
  resultsNames(dds)
  res<-results(dds, pAdjustMethod="BH", contrast=c("hormone1","hormone2"))}