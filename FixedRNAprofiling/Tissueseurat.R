library(Seurat)
#readfilteredmatrix files
tissue.counts <- Read10X(data.dir = "/athena/chenlab/scratch/jjv4001/Tissue/sample_filtered_feature_bc_matrix")
tissue <- CreateSeuratObject(counts = tissue.counts, min.cells = 3, min.features = 200)
metadata<-read.csv("/athena/chenlab/scratch/jjv4001/tissues.csv")
tissue$sample<-metadata$sample
#normalization, scaling and clustering
tissue <- NormalizeData(tissue)
tissue <- FindVariableFeatures(tissue, selection.method = "vst", nfeatures = 2000)
tissue <- ScaleData(tissue, features = rownames(tissue))
tissue <- RunPCA(tissue, features = VariableFeatures(object = tissue))
tissue <- FindNeighbors(tissue, dims = 1:10)
tissue <- FindClusters(tissue, resolution = 0.1)
tissue <- RunUMAP(tissue, dims = 1:10)
Idents(tissue)<-tissue$sample