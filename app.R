---
title: "App.R"
author: "Lay Teng Ang"
date: "11/12/2020"
output: html_document
---

##Installation scClustViz
devtools::install_github("BaderLab/scClustViz")

library(Seurat)
library(rhdf5)
library (patchwork)
library (dplyr)
library(scClustViz)
library(devtools)
library(presto)
install_github('immunogenomics/presto')

h1 <- Read10X("/Users/LTANG/Downloads/h1_filtered_gene_bc_matrix/")
data.h1 <- CreateSeuratObject(h1, project = "h1", min.cells = 3, min.features = 200)
data.h1[["percent.mt"]] <- PercentageFeatureSet(data.h1, pattern = "^MT-")
data.h1 <- subset(data.h1, subset = nFeature_RNA > 3750 & nFeature_RNA < 8000 & percent.mt <14 & percent.mt >7.5)
data.h1 <- NormalizeData(data.h1)
data.h1 <- FindVariableFeatures(data.h1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(data.h1)
data.h1 <- ScaleData(data.h1, features = all.genes)
data.h1 <- RunPCA(data.h1, features = VariableFeatures(object = data.h1))
ElbowPlot(data.h1)

data.h1 <- FindNeighbors(data.h1, dims = 1:10)
data.h1 <- FindClusters(data.h1, resolution = 0.5)
DATA.h1 <- RunUMAP(data.h1, dims = 1:10)
DimPlot(DATA.h1, reduction = "umap")
save(DATA.h1, file = "/Users/LTANG/Downloads/Artery_vein/h1_app.rds")


## Basic Usage----


# if using Seurat, this regex can grab 
# the metadata columns representing cluster results:
your_cluster_columns <- grepl("res[.0-9]+$",
                              names(getMD(DATA.h1)))

your_cluster_results <- getMD(DATA.h1)[,your_cluster_columns]
?CalcAllSCV
sCVdata_list <- CalcAllSCV(
  inD=DATA.h1,
  clusterDF=your_cluster_results,
  assayType='RNA', #specify assay slot of data
  DRforClust="pca",#reduced dimensions for silhouette calc
  exponent=exp(1), #log base of normalized data
  pseudocount=1,
  DRthresh=0.1, #gene filter - minimum detection rate
  testAll=F, #stop testing clusterings when no DE between clusters
  FDRthresh=0.05,
  calcSil=T, #use cluster::silhouette to calc silhouette widths
  calcDEvsRest=T,
  calcDEcombn=T
)

save(DATA.h1,sCVdata_list,file="/Users/LTANG/Downloads/ArteryVein/h1_app.RData")

#shiny----

runShiny(
  filePath="/Users/LTANG/Downloads/ArteryVein/h1_app.RData",
  
  outPath="./",
  # Save any further analysis performed in the app to the
  # working directory rather than library directory.
  
  annotationDB="org.Hs.eg.db",
  # This is an optional argument, but will add annotations.
  
  cellMarkers=list("ESC"=c("MKI67","SOX2","POU5F1",
                                           "NANOG", "KLF4"),
                   "Primitive Streak"=c("TBXT","MIXL1"),
                   "Lateral Mesoderm"="HAND1","ETV1","LMO2","SCL",
                   "Artery Endothelial"="DLL4","CXCR4","ESM1","SOX17",
                   "Pre-Vein Endothelial"="CD144", "CD31",
                   "Vein Endothelial"=c("NR2F2","NT5E"),
                   "Mesenchymal"=c("ACTC1","PDGFRB")
  ),
  # This is a list of canonical marker genes per expected cell type.
  # The app uses this list to automatically annotate clusters.
  
  imageFileType="png"
  #Set the file format of any saved figures from the app.
)

#RSConnect----

install.packages('rsconnect')
library(rsconnect)

SECRET<-'PU0gn+NE3Jt+8G6IdMYQmRS4G7mMV0Uh+DKk0EgL'

rsconnect::setAccountInfo(name='anglab-endothelial',
                          token='D00AA67F92AB6FF8159D36F963DC5379',
                          secret='PU0gn+NE3Jt+8G6IdMYQmRS4G7mMV0Uh+DKk0EgL')

rsconnect::deployApp('/Users/LTANG/Downloads/ArteryVein/')
?deployApp

#Data Package----

devtools::install_github("BaderLab/HumanLiver") 
library(HumanLiver)
viewHumanLiver()

#making own data package
dir.create("/Users/LTANG/Downloads/ArteryVein/inst/",recursive=T)
save(DATA.h1,sCVdata_list,
     file="/Users/LTANG/Downloads/ArteryVein/inst/h1_app.RData")

runShiny("/Users/LTANG/Downloads/ArteryVein/inst/h1_app.RData")




