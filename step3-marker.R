rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(ggsci)
getwd()
setwd('3-cell/')
sce.all=readRDS( "../2-harmony/sce.all_int.rds")
sce.all
Idents(sce.all) <- sce.all@meta.data$RNA_snn_res.0.3
DimPlot(sce.all,label = T)

  celltype=data.frame(ClusterID=0:18,
                      celltype= 'NA' ) 
  #定义细胞亚群     
  celltype[celltype$ClusterID %in% c( 0,2,7,15,16 ),2]='Fibro'
  celltype[celltype$ClusterID %in% c( 4,6 ),2]='Krt'
  celltype[celltype$ClusterID %in% c(12 ),2]='Cyc' 
  
  celltype[celltype$ClusterID %in% c( 5),2]='Myeloid'  
  celltype[celltype$ClusterID %in% c( 3),2]='SMC'  
  celltype[celltype$ClusterID %in% c( 1,9,14,17 ),2]='Endo' 
  celltype[celltype$ClusterID %in% c( 10),2]='LyEndo' 
  celltype[celltype$ClusterID %in% c( 13),2]='Mast'
  celltype[celltype$ClusterID %in% c( 11),2]='Schwan'
  celltype[celltype$ClusterID %in% c( 18),2]='Melan'
  celltype[celltype$ClusterID %in% c( 8),2]='Lymphocyte'
  
  head(celltype)
  celltype
  table(celltype$celltype)
  sce.all@meta.data$celltype = "NA"
  for(i in 1:nrow(celltype)){
    sce.all@meta.data[which(sce.all@meta.data$RNA_snn_res.0.3 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
  table(sce.all@meta.data$celltype)
  sce <- sce.all
  genes_to_check =c("DCN","LUM",
                    "Krt1","Krt5",
                    "Pecam1","Cdh5",
                    "AIF1","LYZ",
                    "PTPRC","CD3D",
                    "CNN1","ACTA2",
                    "S100B","CDH19",
                    "LYVE1","PROX1",
                    "Top2a","Mki67",
                    "CPA3","TPSAB1",
                    "MLANA","DCT"
  )
  
  genes_to_check = lapply(genes_to_check, toupper)
  genes_to_check <- as.character(genes_to_check)
  genes_to_check= unique(genes_to_check)
  th=theme(axis.text.x = element_text(angle = 45, 
                                      vjust = 0.5, hjust=0.5)
           ,text = element_text(size = 24)
           ,axis.text = element_text(size = 24)) 
  Idents(sce) <- sce@meta.data$celltype
  p <- DotPlot(sce, features = unique(genes_to_check),
               assay='RNA'  )+th+
    scale_color_distiller(palette = "Spectral")
  
  p
  
    DimPlot(sce, reduction = "umap", group.by = "celltype",label = F,
            split.by = "group",pt.size = 0.5)+
      theme(text = element_text(size = 20))+ 
      scale_color_brewer(palette = "Paired")
    
    pro = 'cosg_celltype_'
    library(COSG)
    table(Idents(sce))  
    marker_cosg <- cosg(
      sce,
      groups='all',
      assay='RNA',
      slot='data',
      mu=1,
      n_genes_user=100)
    
    
    ## Top10 genes
    library(dplyr)  
    top_10 <- unique(as.character(apply(marker_cosg$names,2,head,10)))
    # width <-0.006*dim(sce)[2];width
    # height <- 0.25*length(top_10)+4.5;height
    
    width <- 15+0.5*length(unique(Idents(sce)));width
    height <- 8+0.1*length(top_10);height
    library(RColorBrewer)
    display.brewer.all()
    brewer.pal(11, 'Paired')
    DoHeatmap( subset(sce,downsample=100), top_10 , 
               size=3,group.colors = c("#B2DF8A","#33A02C","#1F78B4","#CAB2D6","#E31A1C",
                                       "#FFFF99","#6A3D9A","#FB9A99","#A6CEE3","#FDBF6F","#FF7F00")
    )+
      theme(text = element_text(size = 16))+
      scale_fill_gradientn(colors = c("navy","white","firebrick3"))

    
    p_Fib <- FeaturePlot(sce,features = "DCN",cols = c("lightgrey" ,"#DE1F1F"),slot = "data",label.size = 20,pt.size = 1.2) + 
      theme(text = element_text(size = 35),axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
    p_Fib
    p_Endo <- FeaturePlot(sce,features = "PECAM1",cols = c("lightgrey" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 1.2) + 
      theme(text = element_text(size = 35),axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
    p_Endo
    p_Krt <- FeaturePlot(sce,features = "KRT1",cols = c("lightgrey" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 1.2) + 
      theme(text = element_text(size = 35),axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
    p_Krt 
    p_Myeloid <- FeaturePlot(sce,features = "AIF1",cols = c("lightgrey" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 1.2) + 
      theme(text = element_text(size = 35),axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
    p_Myeloid
    p_Lym <- FeaturePlot(sce,features = "PTPRC",cols = c("lightgrey" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 1.2) + 
      theme(text = element_text(size = 35),axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
    p_Lym 
    p_SMC <- FeaturePlot(sce,features = "ACTA2",cols = c("lightgrey" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 1.2) + 
      theme(text = element_text(size = 35),axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
    p_SMC 
    p_Sch <- FeaturePlot(sce,features = "S100B",cols = c("lightgrey" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 1.2) + 
      theme(text = element_text(size = 35),axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
    p_Sch
    p_LyEndo <- FeaturePlot(sce,features = "LYVE1",cols = c("lightgrey" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 1.2) + 
      theme(text = element_text(size = 35),axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
    p_LyEndo
    p_Cyc <- FeaturePlot(sce,features = "MKI67",cols = c("lightgrey" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 1.2) + 
      theme(text = element_text(size = 35),axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
    p_Cyc
    p_Mast <- FeaturePlot(sce,features = "CPA3",cols = c("lightgrey" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 1.2) + 
      theme(text = element_text(size = 35),axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
    p_Mast
    p_Melan <- FeaturePlot(sce,features = "MLANA",cols = c("lightgrey" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 1.2) + 
      theme(text = element_text(size = 35),axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + guides(color=F)
    p_Melan
    
