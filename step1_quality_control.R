rm(list=ls())
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)

###### step1:import data ######   
library(data.table)
dir='../01_cluster/raw data/' 
samples=list.files( dir )
samples 

sceList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro) 
  ct = Read10X(file.path(dir,pro))
  if(length(ct)==2){ct=ct[[1]]}
  sce=CreateSeuratObject( ct, 
                          project = gsub('fewrgewgergre','',pro),
                          min.cells = 5,
                          min.features = 300 ) 
  return(sce)
})
names(sceList)  

samples

sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids = gsub('fewrgewgergre','',samples))

as.data.frame(sce.all@assays$RNA@counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident)

# 
library(stringr)
phe=str_split( sce.all$orig.ident ,'_',simplify = T)
head(phe);dim(phe)

table(phe[,2])
table(phe[,3])


sce.all@meta.data$group=phe[,2]
#sce.all@meta.data$patients=phe[,2]
table(sce.all@meta.data$group)
table(sce.all$group,sce.all$orig.ident)


###### step2:QC ######
dir.create("./1-QC")
setwd("./1-QC")
mito_genes=rownames(sce.all)[grep("^MT-", rownames(sce.all))] 
mito_genes 
sce.all=PercentageFeatureSet(sce.all, "^MT-", col.name = "percent_mito")
fivenum(sce.all@meta.data$percent_mito)

ribo_genes=rownames(sce.all)[grep("^Rp[sl]", rownames(sce.all),ignore.case = T)]
ribo_genes
sce.all=PercentageFeatureSet(sce.all, "^RP[SL]", col.name = "percent_ribo")
fivenum(sce.all@meta.data$percent_ribo)

rownames(sce.all)[grep("^Hb[^(p)]", rownames(sce.all),ignore.case = T)]
sce.all=PercentageFeatureSet(sce.all, "^HB[^(P)]", col.name = "percent_hb")
fivenum(sce.all@meta.data$percent_hb)

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
feats <- c("nFeature_RNA", "nCount_RNA")
p1=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
  NoLegend()
p1
library(ggplot2) 
ggsave(filename="Vlnplot1.pdf",plot=p1)
feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, same.y.lims=T) + 
  scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend()
p2	
ggsave(filename="Vlnplot2.pdf",plot=p2)

p3=FeatureScatter(sce.all, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
ggsave(filename="Scatterplot.pdf",plot=p3)


selected_c <- WhichCells(sce.all, expression = nFeature_RNA > 300)
selected_f <- rownames(sce.all)[Matrix::rowSums(sce.all@assays$RNA@counts > 0 ) > 3]

sce.all.filt <- subset(sce.all, features = selected_f, cells = selected_c)
dim(sce.all) 
dim(sce.all.filt) 

C=sce.all.filt@assays$RNA@counts
dim(C)
C=Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100

C=C[,sample(1:ncol(C),1000)]
most_expressed <- order(apply(C, 1, median), decreasing = T)[50:1]
pdf("TOP50_most_expressed_gene.pdf",width=14)
boxplot(as.matrix(Matrix::t(C[most_expressed, ])),
        cex = 0.1, las = 1, 
        xlab = "% total count per cell", 
        col = (scales::hue_pal())(50)[50:1], 
        horizontal = TRUE)
dev.off()
rm(C)


selected_mito <- WhichCells(sce.all.filt, expression = percent_mito < 20)
selected_ribo <- WhichCells(sce.all.filt, expression = percent_ribo > 3)
selected_hb <- WhichCells(sce.all.filt, expression = percent_hb < 1 )
length(selected_hb)
length(selected_ribo)
length(selected_mito)


sce.all.filt <- subset(sce.all.filt, cells = selected_mito)
sce.all.filt <- subset(sce.all.filt, cells = selected_ribo)
sce.all.filt <- subset(sce.all.filt, cells = selected_hb)
dim(sce.all.filt)

table(sce.all.filt$orig.ident) 


feats <- c("nFeature_RNA", "nCount_RNA")
p1_filtered=VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
  NoLegend()
ggsave(filename="Vlnplot1_filtered.pdf",plot=p1_filtered)

feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2_filtered=VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
  NoLegend()
ggsave(filename="Vlnplot2_filtered.pdf",plot=p2_filtered)


sce.all.filt <- sce.all.filt[!grepl("MALAT1", rownames(sce.all.filt),ignore.case = T), ]

sce.all.filt <- sce.all.filt[!grepl("^MT-", rownames(sce.all.filt),ignore.case = T), ]


dim(sce.all.filt) 

sce.all.filt = NormalizeData(sce.all.filt)
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
sce.all.filt=CellCycleScoring(object = sce.all.filt, 
                              s.features = s.genes, 
                              g2m.features = g2m.genes, 
                              set.ident = TRUE)
p4=VlnPlot(sce.all.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", 
           ncol = 2, pt.size = 0)
ggsave(filename="Vlnplot4_cycle.pdf",plot=p4)

sce.all.filt@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal()
ggsave(filename="cycle_details.pdf" )


dim(sce.all) 

setwd('../')

source('step2-harmony.R')
