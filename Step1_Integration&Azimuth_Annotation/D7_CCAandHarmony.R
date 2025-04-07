#Upload Packages
library(Seurat)
library(dplyr)
library(tidyverse)
library(rhdf5)
library(ggplot2)
library(patchwork)
library(Matrix)
library(harmony)

setwd(dir = "Downloads/CART_GitHub/Data/")  

#Good Data Upload
main_directory_g <- "D7_g/"
D7_Samples_g <- list()
files_g <- list.files(main_directory_g, full.names = TRUE)

#Using only the first 3 h5 files 
files_g <- files_g[1:3]

for (file in files_g) {
  print(paste("Processing file:", file))
  d7_data_g <- Read10X_h5(file)
  
  g1 <- CreateSeuratObject(counts = d7_data_g$`Gene Expression`, project = "Good_d7_p110", min.cells = 3, min.features = 200)
  g1 <- PercentageFeatureSet(g1, pattern = "^MT-", col.name = "percent.mt")
  g1 <- subset(g1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)
  
  g1$id <- basename(file)
  D7_Samples_g[[basename(file)]] <- g1
}


#Haradvhala Data Upload

main_directory_h <-"D7_h/"
folders_h <- list.dirs(path = main_directory_h, full.names = TRUE, recursive = FALSE)
D7_Samples_h <- list()

#Using only the first 3 folders
folders_h <- folders_h[1:3]

for (folder in folders_h) {
  print(paste("Processing folder:", folder))
  data_d7_h <- Read10X(data.dir = folder)
  h1 <- CreateSeuratObject(data_d7_h, min.cells = 3, min.features = 200)
  h1 <- PercentageFeatureSet(h1, pattern = "^MT-", col.name = "percent.mt")
  h1 <- subset(h1, subset = nFeature_RNA > 200 &nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)
  h1$id <- basename(folder)
  D7_Samples_h[[basename(folder)]] <- h1
}

#Louie Data Upload
main_directory_louie <- "D7and8_l/"
folders_louie <- list.dirs(path = main_directory_louie, full.names = TRUE, recursive = FALSE)
D7_Samples_l <- list()

#Using only the first 3 folders
folders_louie <- folders_louie[1:3]

for (folder_l in folders_louie) {
  print(paste("Processing folder:", folder_l))
  
  barcode.path1 <- file.path(folder_l, "barcodes.tsv")
  genes.path1 <- file.path(folder_l, "genes.tsv")
  matrix.path1 <- file.path(folder_l, "matrix.mtx")
  barcodes1 <- readLines(barcode.path1)
  genes1 <- readLines(genes.path1)
  matrix1 <- readMM(file = matrix.path1) %>% as("dgCMatrix")
  rownames(matrix1) <- genes1
  matrix1 <- matrix1[!duplicated(genes1),]
  l1 <- CreateSeuratObject(counts = matrix1)
  l1 <- PercentageFeatureSet(l1, pattern = "^MT-", col.name = "percent.mt")
  l1 <- subset(l1, subset = nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)
  l1$id <- basename(folder_l)
  D7_Samples_l[[basename(folder_l)]] <- l1
}

#Merge Seurat Objects 
x_merge <- D7_Samples_g[[1]]
y_merge <- c(D7_Samples_g[2:3], D7_Samples_h,D7_Samples_l)
obj <- Reduce(function(x, y) merge(x, y), y_merge, init = x_merge)

# Save the object 
write_rds(x = obj, file = "obj_d7.rds")


# Integration

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
ElbowPlot(obj) # The first 10 PCA will be selected

obj <- FindNeighbors(obj)
obj <- FindClusters(obj, cluster.name = "unintegrated_clusters", resolution = 0.7, algorithm = 2)
obj <- RunUMAP(obj, dims = 1:10)

plot1 <- DimPlot(obj, group.by = "id", label = FALSE) + ggtitle("Unintegrated D7") + NoLegend()
plot1

obj_int <- IntegrateLayers(object = obj, method = CCAIntegration, orig.reduction = "pca",
                         new.reduction = "integrated.cca", verbose = FALSE, k.weight = 50)
obj_int <- FindNeighbors(obj_int, reduction = "integrated.cca", dims = 1:10)
obj_int <- FindClusters(obj_int, cluster.name = "cca_clusters", resolution = 0.7, algorithm = 2)
obj_int <- RunUMAP(obj_int, dims = 1:10, reduction = "integrated.cca", reduction.name = "umap.cca")

plot2 <- DimPlot(obj_int,reduction = "umap.cca", group.by = "id") + ggtitle("Integrated CCA D7")+NoLegend()
plot2

#Harmony Integration

obj_har <- NormalizeData(obj)
obj_har <- FindVariableFeatures(obj_har)
obj_har <- ScaleData(obj_har)
obj_har <- RunPCA(obj_har)
ElbowPlot(obj_har) # The first 10 PCA will be selected

obj_har <- RunHarmony(obj_har, group.by.var = "id")

obj_har <- RunUMAP(obj_har, reduction= "harmony", dims = 1:10)
obj_har <- FindNeighbors(obj_har, reduction= "harmony", dims = 1:10)
obj_har <- FindClusters(obj_har)

plot_a_har <- DimPlot(obj_har,group.by = "id") + ggtitle("Harmony samples") 
plot_b_har <- DimPlot(obj_har,group.by = "seurat_clusters", label = TRUE) + ggtitle("Harmony cell clusters") + NoLegend()
plot_a_har + plot_b_har

#### Save Objects
SaveSeuratRds(object = obj, file = "obj_for_azimuth.rds")
SaveSeuratRds(object = obj_int, file = "obj2_for_azimuth.rds")
SaveSeuratRds(object = obj_har, file = "obj_har_for_azimuth.rds")

#### Azimuth Cell Type Annotation

library(Azimuth)

megeredobj = JoinLayers(obj)
megeredobj_2 = JoinLayers(obj_int)
megeredobj_har = JoinLayers(obj_har)


obj_azimuth <- RunAzimuth(megeredobj, reference = "pbmcref")
plot1_azimuth <- DimPlot(obj_azimuth, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend() +ggtitle("unintegrated D7") 

obj_2_azimuth <- RunAzimuth(megeredobj_2, reference = "pbmcref")
plot2_azimuth <- DimPlot(obj_2_azimuth,reduction = "umap.cca", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 , repel = TRUE) + NoLegend()  + ggtitle("CCA D7")

obj_har_azimuth <- RunAzimuth(megeredobj_har, reference = "pbmcref")
plot3_azimuth <- DimPlot(obj_har_azimuth, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 , repel = TRUE) + NoLegend() + ggtitle("Harmony D7")

plot1_azimuth + plot2_azimuth + plot3_azimuth

#### SAVE Azimuth annotated RDS

SaveSeuratRds(object = obj_azimuth, file = "azimuth_obj_D7_unintegrated.rds")
SaveSeuratRds(object = obj_2_azimuth, file = "azimuth_obj_D7_CCA.rds")
SaveSeuratRds(object = obj_har_azimuth, file = "azimuth_obj_D7_Harmony.rds")

