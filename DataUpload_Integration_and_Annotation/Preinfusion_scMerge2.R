library(SingleCellExperiment)
library(scMerge)
library(Matrix)
library(scuttle)
library(Seurat)

setwd(dir = "Downloads/CART_GitHub/Data/")  

dirs <- list(
  Louie = "Preinfusion_l",
  Haradhvala = "Preinfusion_h"
)


process_directory <- function(parent_dir, id_prefix) {
  patient_folders <- list.dirs(parent_dir, recursive = FALSE)
  sce_sublist <- list()
  
  for (folder in patient_folders[1:3]) {
    barcodes_file <- list.files(folder, pattern = "barcodes", full.names = TRUE)
    features_file <- list.files(folder, pattern = "features|genes", full.names = TRUE)
    matrix_file <- list.files(folder, pattern = "matrix", full.names = TRUE)
    
    barcodes <- read.delim(barcodes_file, header = FALSE)
    features <- read.delim(features_file, header = FALSE)
    matrix <- as(readMM(matrix_file), "dgCMatrix")  # Convert to dgCMatrix
    
    if (id_prefix == "Louie") {
      rownames(matrix) <- features$V1
    } else if (id_prefix == "Haradhvala") {
      rownames(matrix) <- features$V2
    }
    colnames(matrix) <- barcodes$V1
    
    # Create SingleCellExperiment object
    sce <- SingleCellExperiment(
      assays = list(counts = matrix),
      rowData = data.frame(gene = features$V1),
      colData = data.frame(cell = barcodes$V1)
    )
    
    # Filter based on quality control metrics
    sce <- sce[rowSums(counts(sce) > 0) > 3, ]  # Filter genes detected in fewer than 3 cells
    sce <- sce[, colSums(counts(sce)) > 200 & colSums(counts(sce)) < 4000]  # Filter cells by feature counts
    mitochondrial_genes <- grepl("^MT-", rowData(sce)$gene)
    sce$percent.mt <- colSums(counts(sce)[mitochondrial_genes, ]) / colSums(counts(sce)) * 100
    sce <- sce[, sce$percent.mt < 10]  # Filter cells with high mitochondrial content
    
    # Add batch ID and metadata
    sce$sample <- paste0(id_prefix, "_", basename(folder))
    colData(sce)$sample_type <- rep("Preinfusion", ncol(sce))
    colData(sce)$condition <- rep("Baseline", ncol(sce))
    
    sce_sublist[[basename(folder)]] <- sce
  }
  
  return(sce_sublist)
}

# Function to process each directory and create SingleCellExperiment objects

sce_list <- list()

for (dir_name in names(dirs)) {
  sce_list[[dir_name]] <- process_directory(dirs[[dir_name]], id_prefix = dir_name)
}

## Combine all SCE objects into a single list
sce_list <- do.call(c, sce_list)
#Labels
batch_labels <- names(sce_list) 

# Update the `gene` column in `rowData` for all SCE objects
for (name in names(sce_list)) {
  sce <- sce_list[[name]]
  
  # Set the `gene` column in `rowData` to match the row names
  rowData(sce)$gene <- rownames(sce)
  
  sce_list[[name]] <- sce
}

common_genes <- Reduce(intersect, lapply(sce_list, rownames))
sce_list <- lapply(sce_list, function(sce) sce[common_genes, ])
combined_sce <- do.call(cbind, sce_list)
batch_labels <- rep(names(sce_list), times = sapply(sce_list, ncol))

# Add batch information to colData of the combined object
colData(combined_sce)$batch <- batch_labels

logcounts(combined_sce) <- log1p(counts(combined_sce))

colnames(combined_sce@assays@data$logcounts) <- make.unique(colnames(combined_sce@assays@data$logcounts))

combined_sce <- scater::runUMAP(combined_sce, ncomponents = 2)
plot1 <- scater::plotUMAP(combined_sce, colour_by = "sample")+
  ggplot2::ggtitle("scMerge2 Integrated Cells")

plot1

# Save the UMAP plot if you want
ggplot2::ggsave("scMerge2_Unintegrated_umap.png", plot = plot1, width = 8, height = 6, dpi = 300)

# Perform integration using scMerge2
integrated_sce <- scMerge2(
  exprsMat = combined_sce@assays@data$logcounts,
  batch = colData(combined_sce)$batch
)

#add scMerge2 to the sce object
assay(combined_sce, "scMerge2", withDimnames = FALSE) <- integrated_sce$newY

#   combined_sce = scater::runPCA(combined_sce, exprs_values = 'scMerge2')
#   scater::plotPCA(combined_sce, colour_by = 'sample')

combined_sce <- scater::runUMAP(combined_sce, exprs_values = "scMerge2", ncomponents = 2)
plot2 <- scater::plotUMAP(combined_sce, colour_by = "sample")+
  ggplot2::ggtitle("scMerge2 Integrated Cells")

plot2

# Save the UMAP plot if you want
ggplot2::ggsave("scMerge2_Integrated_umap.png", plot = plot2, width = 8, height = 6, dpi = 300)

#save in current directory
saveRDS(combined_sce, file = "scMerge2_Pre.rds")



# Azimuth Annotation
library(Azimuth)

# Check for duplicated cell names
table(duplicated(colnames(combined_sce)))

# If duplicates exist, forcibly make them unique
colnames(combined_sce) <- make.unique(colnames(combined_sce))

# Now try converting again
scMerge_seurat <- as.Seurat(combined_sce)

#scMerge_seurat_merged <- JoinLayers(scMerge_seurat)

obj_azimuth <- RunAzimuth(scMerge_seurat, assay = "originalexp", reference = "pbmcref")
obj_azimuth

plot3 <- DimPlot(obj_azimuth, reduction = "UMAP", group.by = "sample", raster = FALSE)+NoLegend()
plot4 <- DimPlot(obj_azimuth, reduction = "UMAP", group.by = "predicted.celltype.l1", raster = FALSE)

plot3 + plot4

# But first add the missing metadata columns that will be needed for the analysis
metadata <- obj_azimuth@meta.data

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
obj_azimuth@meta.data <- metadata

