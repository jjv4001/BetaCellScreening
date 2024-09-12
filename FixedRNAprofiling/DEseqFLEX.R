library(tidyverse)
library(DESeq2)
library(Seurat)
tissues<-readRDS("/athena/chenlab/scratch/jjv4001/tissuesnew.RDS")
Idents(tissues)<-tissues$sample
FemaleTreated<-subset(tissues, sample=="Female Treated")
MaleTreated<-subset(tissues, sample=="Male Treated")
num_cells_FemaleTreated <- ncol(FemaleTreated)
FemaleTreated$cell_label <- paste0("cell", 1:num_cells_FemaleTreated)
num_cells_MaleTreated <- ncol(MaleTreated)
MaleTreated$cell_label <- paste0("cell", 1:num_cells_MaleTreated)
combined<-merge(FemaleTreated,MaleTreated)
counts<-AggregateExpression(combined, group.by= c("cell_label","sample"),slot= "counts", return.seurat = FALSE)
counts<-counts$RNA
#order columns in counts matrix according to names 
suffix <- gsub(".*_", "", colnames(counts))
sorted_colnames <- colnames(counts)[order(suffix)]
counts <- counts[, sorted_colnames]
#create colData corresponding to each cell and hormone to be distinguished (col1=cell1_hormone1, col2=hormone1 or col1=cell100_hormone2, col2=hormone2)
colData<-data.frame(samples=colnames(counts))
rownames(colData) <- NULL
colData$sample<-gsub('.*_','',colData$samples)
colData <- colData %>% 
  column_to_rownames(var = 'samples')
dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = colData, design = ~sample)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds<-estimateSizeFactors(dds, type="poscount")
dds<-DESeq(dds)
resultsNames(dds)
res<-results(dds, pAdjustMethod="BH", contrast=c("sample","Female Treated","Male Treated"))
write.csv(res,"/athena/chenlab/scratch/jjv4001/FemalevsMaleTreated.csv")
