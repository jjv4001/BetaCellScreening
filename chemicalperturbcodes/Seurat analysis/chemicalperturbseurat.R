library(Seurat)
#readfilteredmatrix files
endoc1.counts <- Read10X(data.dir = "/athena/chenlab/scratch/jjv4001/source/Jeya_1/sample_filtered_feature_bc_matrix")
endoc1 <- CreateSeuratObject(counts = endoc1.counts, min.cells = 3, min.features = 200)
endoc2.counts <- Read10X(data.dir = "/athena/chenlab/scratch/jjv4001/source/Jeya_2/sample_filtered_feature_bc_matrix")
endoc2 <- CreateSeuratObject(counts = endoc2.counts, min.cells = 3, min.features = 200)
endoc3.counts <- Read10X(data.dir = "/athena/chenlab/scratch/jjv4001/source/Jeya_3/sample_filtered_feature_bc_matrix")
endoc3 <- CreateSeuratObject(counts = endoc3.counts, min.cells = 3, min.features = 200)
endoc<-merge(x=endoc1, y=c(endoc2, endoc3), add.cell.ids=c("endoc1","endoc2","endoc3"))
metadata<-read.csv("/athena/chenlab/scratch/jjv4001/endocchemperturbmetadata.csv")
endoc$hormone<-metadata$hormone
#normalization, scaling and clustering
endoc <- NormalizeData(endoc)
endoc <- FindVariableFeatures(endoc, selection.method = "vst", nfeatures = 2000)
endoc <- ScaleData(endoc, features = rownames(endoc))
endoc <- RunPCA(endoc, features = VariableFeatures(object = endoc))
endoc <- FindNeighbors(endoc, dims = 1:10)
endoc <- FindClusters(endoc, resolution = 0.18)
endoc <- RunUMAP(endoc, dims = 1:10)
Idents(endoc)<-endoc$hormone