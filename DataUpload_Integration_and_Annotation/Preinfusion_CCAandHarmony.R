#Aim: is to upload Data For Pre-infusion Total 24 Sample
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)

#Set thw orking directory where you have your data
setwd(dir = "Downloads/CART_GitHub/Data/")  

#Haradvhala Data Upload 
main_directory <- "Preinfusion_h/"
folders_h <- list.dirs(path = main_directory, full.names = TRUE, recursive = FALSE)
Preinf_Samples_h <- list()

#Selecting only 3 samples
folders_h <- folders_h[1:3]

# Seurat Object Creation for Haradvhala
for (folder in folders_h) {
  print(paste("Processing folder:", folder))
  data_preinf_h <- Read10X(data.dir = folder) #Each folder contains barcodes, genes and matrix files
  
  h_pre <- CreateSeuratObject(data_preinf_h, min.cells = 3, min.features = 200) #Creating Seurat Object
  h_pre <- PercentageFeatureSet(h_pre, pattern = "^MT-", col.name = "percent.mt") #Calculating percentage of mitochondrial genes
  h_pre <- subset(h_pre, subset = nFeature_RNA > 200 &nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10) #Filtering cells based on number of genes, counts and percentage of mitochondrial genes
  h_pre$sample <- basename(folder)
  Preinf_Samples_h[[basename(folder)]] <- h_pre
}
print(length(Preinf_Samples_h))


#Louie Data Upload
main_directory_louie <- "Preinfusion_l/"
folders_louie <- list.dirs(path = main_directory_louie, full.names = TRUE, recursive = FALSE)

#Only Three Samples
folders_louie <- folders_louie[1:3]
Preinf_Samples_l <- list()

# Seurat Object Creation for Louie
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
  l1 <- subset(l1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)
  l1$sample <- basename(folder_l)
  Preinf_Samples_l[[basename(folder_l)]] <- l1
}

print(length(Preinf_Samples_l))

#Merge Seurat Objects into One
x_merge <- Preinf_Samples_h[[1]]
y_merge <- c(Preinf_Samples_h[2:3], Preinf_Samples_l)
obj <- Reduce(function(x, y) merge(x, y), y_merge, init = x_merge)

View(obj@meta.data) #You can investigate the metadata of the object if needed

# But first add the missing metadata columns that will be needed for the analysis
metadata <- obj@meta.data

metadata$response <- NA
metadata$Sex <- NA
metadata$Age <- NA

unique_ids <- unique(metadata$sample)
unique_ids

metadata$response[metadata$sample == "GSE224550_P1_PBMC_Preinfusion"] <- "CR"
metadata$response[metadata$sample == "GSE224550_P2_PBMC_Preinfusion"] <- "PD"
metadata$response[metadata$sample == "GSE224550_P7_PBMC_Preinfusion"] <- "CR"
metadata$response[metadata$sample == "GSE224550_P8_PBMC_Preinfusion"] <- "CR"

#Haradhvala  ###############################################

metadata$response[metadata$sample == "Patient10-Baseline"] <- "PD"
metadata$response[metadata$sample == "Patient11-Baseline"] <- "CR"
metadata$response[metadata$sample == "Patient12-Baseline"] <- "CR"
metadata$response[metadata$sample == "Patient13-Baseline"] <- "CR"
metadata$response[metadata$sample == "Patient14-Baseline"] <- "PD"
metadata$response[metadata$sample == "Patient15-Baseline"] <- "CR"
metadata$response[metadata$sample == "Patient17-Baseline"] <- "CR"
metadata$response[metadata$sample == "Patient18-Baseline"] <- "PD"
metadata$response[metadata$sample == "Patient19-Baseline"] <- "CR"
metadata$response[metadata$sample == "Patient20-Baseline"] <- "PD"
metadata$response[metadata$sample == "Patient21-Baseline"] <- "CR"
metadata$response[metadata$sample == "Patient22-Baseline"] <- "CR"
metadata$response[metadata$sample == "Patient23-Baseline"] <- "PD"
metadata$response[metadata$sample == "Patient24-Baseline"] <- "PD"
metadata$response[metadata$sample == "Patient25-Baseline"] <- "PD"
metadata$response[metadata$sample == "Patient30-Baseline"] <- "CR"
metadata$response[metadata$sample == "Patient31-Baseline"] <- "PD"
metadata$response[metadata$sample == "Patient6-Baseline"] <- "PD"
metadata$response[metadata$sample == "Patient8-Baseline"] <- "CR"
metadata$response[metadata$sample == "Patient9-Baseline"] <- "PD"


# Update 'Sex' based on sample
#Louie
metadata$Sex[metadata$sample == "GSE224550_P1_PBMC_Preinfusion"] <- "F"
metadata$Sex[metadata$sample == "GSE224550_P2_PBMC_Preinfusion"] <- "M"
metadata$Sex[metadata$sample == "GSE224550_P7_PBMC_Preinfusion"] <- "M"
metadata$Sex[metadata$sample == "GSE224550_P8_PBMC_Preinfusion"] <- "M"

# Haradhvala
metadata$Sex[metadata$sample == "Patient10-Baseline"] <- "F"
metadata$Sex[metadata$sample == "Patient11-Baseline"] <- "M"
metadata$Sex[metadata$sample == "Patient12-Baseline"] <- "M"
metadata$Sex[metadata$sample == "Patient13-Baseline"] <- "M"
metadata$Sex[metadata$sample == "Patient14-Baseline"] <- "M"
metadata$Sex[metadata$sample == "Patient15-Baseline"] <- "M"
metadata$Sex[metadata$sample == "Patient17-Baseline"] <- "M"
metadata$Sex[metadata$sample == "Patient18-Baseline"] <- "M"
metadata$Sex[metadata$sample == "Patient19-Baseline"] <- "M"
metadata$Sex[metadata$sample == "Patient20-Baseline"] <- "F"
metadata$Sex[metadata$sample == "Patient21-Baseline"] <- "M"
metadata$Sex[metadata$sample == "Patient22-Baseline"] <- "M"
metadata$Sex[metadata$sample == "Patient23-Baseline"] <- "M"
metadata$Sex[metadata$sample == "Patient24-Baseline"] <- "M"
metadata$Sex[metadata$sample == "Patient25-Baseline"] <- "F"
metadata$Sex[metadata$sample == "Patient30-Baseline"] <- "F"
metadata$Sex[metadata$sample == "Patient31-Baseline"] <- "F"
metadata$Sex[metadata$sample == "Patient6-Baseline"] <- "M"
metadata$Sex[metadata$sample == "Patient8-Baseline"] <- "M"
metadata$Sex[metadata$sample == "Patient9-Baseline"] <- "F"

# Update 'Age' based on sample
# Louie
metadata$Age[metadata$sample == "GSE224550_P1_PBMC_Preinfusion"] <- "25"
metadata$Age[metadata$sample == "GSE224550_P2_PBMC_Preinfusion"] <- "66"
metadata$Age[metadata$sample == "GSE224550_P7_PBMC_Preinfusion"] <- "47"
metadata$Age[metadata$sample == "GSE224550_P8_PBMC_Preinfusion"] <- "30"

# Haradhvala
metadata$Age[metadata$sample == "Patient10-Baseline"] <- "66"
metadata$Age[metadata$sample == "Patient11-Baseline"] <- "57"
metadata$Age[metadata$sample == "Patient12-Baseline"] <- "55"
metadata$Age[metadata$sample == "Patient13-Baseline"] <- "70"
metadata$Age[metadata$sample == "Patient14-Baseline"] <- "62"
metadata$Age[metadata$sample == "Patient15-Baseline"] <- "38"
metadata$Age[metadata$sample == "Patient17-Baseline"] <- "48"
metadata$Age[metadata$sample == "Patient18-Baseline"] <- "71"
metadata$Age[metadata$sample == "Patient19-Baseline"] <- "67"
metadata$Age[metadata$sample == "Patient20-Baseline"] <- "79"
metadata$Age[metadata$sample == "Patient21-Baseline"] <- "52"
metadata$Age[metadata$sample == "Patient22-Baseline"] <- "61"
metadata$Age[metadata$sample == "Patient23-Baseline"] <- "70"
metadata$Age[metadata$sample == "Patient24-Baseline"] <- "77"
metadata$Age[metadata$sample == "Patient25-Baseline"] <- "64"
metadata$Age[metadata$sample == "Patient30-Baseline"] <- "34"
metadata$Age[metadata$sample == "Patient31-Baseline"] <- "76"
metadata$Age[metadata$sample == "Patient6-Baseline"] <- "77"
metadata$Age[metadata$sample == "Patient8-Baseline"] <- "65"
metadata$Age[metadata$sample == "Patient9-Baseline"] <- "56"

#Update the metadata
obj@meta.data <- metadata

# Seurat Data Workflow
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
ElbowPlot(obj) # The first 15 PCA will be selected

obj <- FindNeighbors(obj)
#In FindClusters function resolution is a parameter to control the number of clusters, while algorithm is the clustering algorithm to be used
obj <- FindClusters(obj, cluster.name = "unintegrated_pre", resolution = 0.7, algorithm = 2) 
obj <- RunUMAP(obj, dims = 1:15)
#Plotting the UMAP
plot1 <- DimPlot(obj, group.by = "sample") + ggtitle("unintegrated_preinfusion (Samples)")
plot1

#Integration of the samples using Seruat's Canonical Correlation Analysis (CCA) method

obj_cca <- IntegrateLayers(object = obj, method = CCAIntegration, orig.reduction = "pca",
                         new.reduction= "integrated.cca", verbose = FALSE)
ElbowPlot(obj_cca) # The first 15 PCA will be selected

obj_cca <- FindNeighbors(obj_cca, reduction = "integrated.cca", dims = 1:15)
obj_cca <- FindClusters(obj_cca, cluster.name = "cca_clusters", resolution = 0.7, algorithm = 2) #resolution is a parameter to control the number of clusters, while algorthim is the clustering algorithm to be used
obj_cca <- RunUMAP(obj_cca, dims = 1:15, reduction = "integrated.cca", reduction.name = "umap.cca")
plot2 <- DimPlot(obj_cca,reduction = "umap.cca", group.by = "sample") + ggtitle("cca_preinfusion")
plot2

#Harmony Integration
library(harmony)

obj_har <- NormalizeData(obj)
obj_har <- FindVariableFeatures(obj_har)
obj_har <- ScaleData(obj_har)
obj_har <- RunPCA(obj_har)
ElbowPlot(obj_har) # The first 15 PCA will be selected

obj_har <- RunHarmony(obj_har, group.by.vars = "sample", dims.use = 1:15)
obj_har <- FindNeighbors(obj_har, reduction= "harmony", dims = 1:15)
obj_har <- FindClusters(obj_har, resolution = 0.7, algorithm = 2) #resolution is a parameter to control the number of clusters, while algorthim is the clustering algorithm to be used
obj_har <- RunUMAP(obj_har, reduction= "harmony", dims = 1:15)

# Plotting the UMAPs
plot1_har <- DimPlot(obj_har,group.by = "sample", label = TRUE, repel = TRUE) + ggtitle("Harmony (Samples)") + NoLegend()
plot2_har <- DimPlot(obj_har,group.by = "seurat_clusters", label = TRUE) + ggtitle("Harmony (Cells)") + NoLegend()
plot1_har
plot2_har

#### Save the Seurat Objects for Later if needed ###
SaveSeuratRds(object = obj, file = "obj_un_4_azimuth.rds")
SaveSeuratRds(object = obj_cca, file = "obj_cca_4_azimuth.rds")
SaveSeuratRds(object = obj_har, file = "obj_har_4_azimuth.rds")

#### Cell Type Annotation using Azimuth tool
library(Azimuth)

#obj <- readRDS("obj_un_4_azimuth.rds")
#obj_cca <- readRDS("obj_cca_4_azimuth.rds")
#obj_har <- readRDS("obj_har_4_azimuth.rds")

#Merging the layers of the Seurat object
megeredobj = JoinLayers(obj)
megeredobj_2 = JoinLayers(obj_cca)
megeredobj_har = JoinLayers(obj_har)

#Running Azimuth for cell type annotation
obj_azimuth <- RunAzimuth(megeredobj, reference = "pbmcref")
obj_cca_azimuth <- RunAzimuth(megeredobj_2, reference = "pbmcref")
obj_har_azimuth <- RunAzimuth(megeredobj_har, reference = "pbmcref")


#Plotting the UMAPs with cell type annotations
plot1_azimuth <- DimPlot(obj_azimuth, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend() +ggtitle("unintegrated Preinfusion") 
plot2_azimuth <- DimPlot(obj_cca_azimuth,reduction = "umap.cca", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 , repel = TRUE) + NoLegend()  + ggtitle("CCA Preinfusion")
plot3_azimuth <- DimPlot(obj_har_azimuth, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 , repel = TRUE) + NoLegend() + ggtitle("Harmony Preinfusion")

plot1_azimuth + plot2_azimuth + plot3_azimuth


#### SAVE Final Seurat objects with Azimuth Annotations ####
SaveSeuratRds(object = obj_azimuth, file = "azimuth_obj_Preinfusion_unintegrated.rds")
SaveSeuratRds(object = obj_2_azimuth, file = "azimuth_obj_Preinfusion_CCA.rds")
SaveSeuratRds(object = obj_har_azimuth, file = "azimuth_obj_Preinfusion_Harmony.rds")


