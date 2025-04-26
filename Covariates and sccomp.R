#Covariates and sccomp

#libraries import
library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(forcats)
library(sccomp)
library(cmdstanr)

#upload the rds Seurat Object
cca_pre <- readRDS("/Users/mohammedkhattab/Downloads/CART_GitHub/Data/obj_har_for_azimuth.rds")

#Create a list of floats corresponding different clustering resolutions 
resolutions <- seq(0.1, 1.2, by = 0.4)


#Create a loop in which each time a new clustering column 
#   will be generated each with a different resolution, each 
#   will be saved as new rds file the number of the created clusters.

for (res in resolutions) {  
  cluster_result <- FindClusters(object = cca_pre, resolution = res, algorithm = 2)

  result_name <- paste("cca_pre_", gsub("\\.", "", res), sep = "")
  assign(result_name, cluster_result)
  
  saveRDS(object = cluster_result, file = paste(result_name, ".rds", sep = ""))
}


#Another way that may sound more convenient is to only
#   have a single file containing the clustering columns.

for (res in resolutions) {
  res_label <- gsub("\\.", "", as.character(res))
  cluster_col_name <- paste0("clusters_res_", res_label)

  cca_pre <- FindClusters(
    object = cca_pre,
    resolution = res,
    algorithm = 2,
    verbose = FALSE
  )

}


# Save the updated Seurat object with all clustering columns
saveRDS(cca_pre, file = "Preinfusion_cca_clusters.rds")



#Before jumping to sccomp a step we found interesting and equally important is the clustering of Azimuth,
#since Azimuth creates clusters of cells based on their expresssion, why not explore these also.

sccomp_result_l1 <- cca_pre |>
  sccomp_estimate(
    formula_composition = ~ percent.mt,
    .sample = id,
    .cell_group = predicted.celltype.l1, #Different cell groups each sccomp run
    cores = 1, verbose = TRUE , variational_inference = TRUE
  ) |>
  sccomp_remove_outliers(cores = 1) |>
  sccomp_test()

# Save the sccomp result as a boxplot image
boxplot_l1 <- sccomp_boxplot(sccomp_result_l1, factor = "response")
png(filename = "cca_predicted.celltype.l1.png", width = 7, height = 7, units = "in", res = 300)
print(boxplot_l1)
dev.off()


## In the second step following this would be the  statstical analysis using sccomp

#Load the sccomp library if each object has a different clustering resolution column
seurat_list <- list(
  cca_09 = readRDS("cca_pre_09.rds"), 
  cca_05 = readRDS("cca_pre_05.rds"),
  cca_01 = readRDS("cca_pre_01.rds")
)

for (obj_name in names(seurat_list)) {
  obj <- seurat_list[[obj_name]]

  cat("Processing:", obj_name, "\n")
  metadata <- obj@meta.data
  metadata <- metadata %>%
    mutate(
      percent.mt = as.numeric(as.character(percent.mt)),
      RNA_snn_res.0.9 = as.character(RNA_snn_res.0.9),
      id = as.character(id),                              #sccomp required columns to be in character
      seurat_clusters = as.character(seurat_clusters),
      Sex = as.character(Sex),
      Age = as.numeric(as.character(Age)) #and like this in the case of a continuous variable
  )

  obj@meta.data <- metadata

  sccomp_result <- obj |>
    sccomp_estimate(
      formula_composition = ~ response + Age + Sex, #determine the covariates to be used
      .sample = id, 
      .cell_group = seurat_clusters, 
      cores = 1, verbose = TRUE
    ) |>
    sccomp_remove_outliers(cores = 1) |>
    sccomp_test()

  sccomp_result_simplified <- sccomp_result %>%
    select(-count_data) #Remove the count data column to simplify the table

  #Save the table for further analysis
  write.csv(sccomp_result_simplified, file = paste0(obj_name, "_response.csv"), row.names = FALSE)

  boxplot <- sccomp_boxplot(sccomp_result, factor = "response") +
    ggtitle("CCA Pre response Sex and Age")

  png(filename = paste0("boxplot_responseSexAge", obj_name,".png"), width = 7, height = 7, units = "in", res = 300)
  print(boxplot)
  dev.off()

}

