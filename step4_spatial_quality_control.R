rm(list = ls())
library("Seurat")
library("ggplot2")
library("patchwork")
library("dplyr")
library("cowplot")

B1 <- Load10X_Spatial("raw data/b1/Spatial+h5/",slice = "B1")
B1@meta.data$orig.ident <- "B1"
B2 <- Load10X_Spatial("raw data/b2/Spatial+h5/",slice = "B2")
B2@meta.data$orig.ident <- "B2"
C1 <- Load10X_Spatial("raw data/c1/Spatial+h5/",slice = "C1")
C1@meta.data$orig.ident <- "C1"
C2 <- Load10X_Spatial("raw data/c2/Spatial+h5/",slice = "C2")
C2@meta.data$orig.ident <- "C2"


B1 <- SCTransform(B1, assay = "Spatial", return.only.var.genes = FALSE, verbose = T)
B2 <- SCTransform(B2, assay = "Spatial", return.only.var.genes = FALSE, verbose = T)
C1 <- SCTransform(C1, assay = "Spatial", return.only.var.genes = FALSE, verbose = T)
C2 <- SCTransform(C2, assay = "Spatial", return.only.var.genes = FALSE, verbose = T)

samples <- c("B1", "B2", "C1","C2")
spatial <- merge(x = B1,y = c(B2,C1,C2),add.cell.ids = samples,project = "keloid")
DefaultAssay(spatial)
VariableFeatures(spatial) <- c(VariableFeatures(B1),VariableFeatures(B2),VariableFeatures(C1),VariableFeatures(C2))
spatial <- RunPCA(spatial, assay = "SCT", verbose = T)
spatial <- RunUMAP(spatial, reduction = "pca", dims = 1:15)
spatial <- FindNeighbors(spatial, reduction = "pca", dims = 1:15)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  spatial=FindClusters(spatial, 
                      resolution = res, algorithm = 1)
}
colnames(spatial@meta.data)
apply(spatial@meta.data[,grep("SCT_snn",colnames(sce.all@meta.data))],2,table)
spatial <- FindClusters(spatial)

SpatialDimPlot(spatial,label = F,group.by = "SCT_snn_res.0.2",
               images = c("B1","C1","B2","C2"),
               pt.size.factor = 1.7)+scale_fill_futurama()
