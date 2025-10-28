#### Clean code for revision and open source code
## Vincent Gélinas
## Émeric Texeraud
## Lou Nielly Thibault
## Charles Joly Beauparlant

## October 14th 2025

library(org.Hs.eg.db)
library(AnnotationDbi)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(ComplexHeatmap)

## Repository -- variable 
# Change these variables as needed to fit your own architecture (these should be created before hand)
# Mulitple figures should be generated with this code, depending on which sections you are using. You can split the different figures on different path depending on the theme

path = "./"
plot_path = "plots/"
results_path = "results/"
data_path = "data/"
our_data = "raw_data/"

### Reading the data and first analysis on Smajic, Kamath, and Feleke #### 

all_data <- readRDS(paste0(our_data, "all_data_raw.rds"))
feleke <- all_data[["feleke"]]
saveRDS(feleke, paste0(our_data, "feleke_raw.rds"))
rm(feleke)
kamath <- all_data[["kamath"]]
saveRDS(kamath, paste0(our_data, "kamath_raw.rds"))
rm(kamath)
smajic <- all_data[["smajic"]]
saveRDS(smajic, paste0(our_data, "smajic_raw.rds"))
rm(smajic)
rm(all_data)
gc()

for (i in c("smajic", "feleke", "kamath")) {
scdata <- readRDS(paste0(our_data, i, "_raw.rds"))
   genes_count <- rowSums(scdata@assays[["RNA"]]@counts > 0)
   genes_filter <- names(genes_count[genes_count >= 10])
   scdata <- subset(scdata, subset = (nFeature_RNA >= 1000))
   scdata <- subset(scdata, features = genes_filter)

   scdata <- Seurat::NormalizeData(scdata)
   scdata <- Seurat::FindVariableFeatures(scdata)
   scdata <- Seurat::ScaleData(scdata)
   scdata <- Seurat::RunPCA(scdata)
   scdata <- Seurat::RunUMAP(scdata, reduction = "pca", dims = 1:30)
   scdata <- Seurat::FindNeighbors(scdata, reduction = "pca", dims = 1:30)
   scdata <- Seurat::FindClusters(scdata, resolution = 0.2)

   scdata <- Seurat::DietSeurat(scdata, assays = "RNA", dimreducs = c("umap", "pca"), counts = TRUE, data = TRUE, scale.data = FALSE)

   saveRDS(scdata, paste0(our_data, i, "_preprocessed.rds"))
   rm(scdata)
   gc()
}

### The Seurat objects are now worked, and the base clusterisation was performed. 

### Working Yang et al., ####

library(Matrix)

control_list <- c("GSM4848456_c15_12", # Choroid plexus non-viral
                  "GSM4848457_c15_17", # Choroid plexus non-viral
                  "GSM4848458_c16_18", # Choroid plexus non-viral
                  "GSM4848459_c16_23", # Choroid plexus non-viral
                  "GSM4848460_c16_24", # Choroid plexus non-viral
                  "GSM4848461_c16_25") # Choroid plexus non-viral

seurat_list <- list()


for (sample in control_list) {
  
  matrix_file <- paste0(our_data, "yang/", sample, "_matrix.mtx")
  barcodes_file <- paste0(our_data, "yang/", sample, "_barcodes.tsv")
  features_file <- paste0(our_data, "yang/", sample, "_features.tsv")
  
 
  expr_matrix <- readMM(matrix_file)  
  barcodes <- read.delim(barcodes_file, header = FALSE)  
  features <- read.delim(features_file, header = FALSE)  
  
  
  colnames(expr_matrix) <- barcodes$V1
  rownames(expr_matrix) <- features$V1
  
  
  seurat_obj <- CreateSeuratObject(counts = expr_matrix, project = sample)
  
  seurat_obj[["orig.ident"]] <- sample
  
  seurat_list[[sample]] <- seurat_obj
}

completed_seurat <- merge(seurat_list[[1]], y = c(seurat_list[[2]],
                                                  seurat_list[[3]],
                                                  seurat_list[[4]],
                                                  seurat_list[[5]],
                                                  seurat_list[[6]]), 
                          add.cell.ids = control_list, project = "Yang")

saveRDS(object = completed_seurat, file = paste0(our_data, "plexus_choroird_control_yang.rds"))

yang <- readRDS(paste0(our_data, "plexus_choroird_control_yang.rds"))

## In-house packages. This function wraps all the important Seurat step for clusterisation. See methodology for better informations.
SingleCell::analyze_integrated(yang, organism = "human", assay = "RNA", force_report = TRUE, file_name = "yang_control_samples.rds", force_DE = FALSE,
                               save_path = paste0(our_data, "analyze_integrated/"))

### Quality control ####

column <- c("nCount_RNA", "nFeature_RNA")

### Smajic

for (analyzed_column in column) {
  df <- data.frame(
    n = smajic@meta.data[[analyzed_column]]
  )
  
  table <- summary(smajic@meta.data[[analyzed_column]])
  
  table <- data.frame(
    Statistic = names(table),
    Value = as.numeric(table)
  )
  
  title_path <- paste0("smajic et al., worked cells; n = ", length(df$n))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = n)) + 
    ggplot2::geom_histogram(bins = 40, color = "grey50", fill = "grey78") + 
    ggplot2::labs(title = title_path, 
                  x = analyzed_column, 
                  y = "count") + 
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[4])), color = "black", linetype = "dashed", size = 0.4) + ## Mean
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[3])), color = "#5FCAE3", linetype = "dashed", size = 0.4) +  ## Median
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[2])), color = "#FF9EBC", linetype = "dashed", size = 0.4) + ## 1st Quant
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[5])), color = "#FF9EBC", linetype = "dashed", size = 0.4) +  ## 3rd Quant
    theme_bw()
  
  graph_name <- paste("smajic", analyzed_column,  "summary.png", sep = "_")
  ggplot2::ggsave(graph_name, plot = p,
                  device = "png", path = plot_path, dpi = 200, width = 1600,
                  height = 700, units = "px")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = n)) + 
    ggplot2::geom_histogram(bins = 40, color = "grey50", fill = "grey78") + 
    ggplot2::labs(title = title_path, 
                  x = analyzed_column, 
                  y = "count") + 
    theme_bw()
  
  graph_name <- paste("smajic", analyzed_column,  "summary_without_stats.png", sep = "_")
  ggplot2::ggsave(graph_name, plot = p,
                  device = "png", path = plot_path, dpi = 200, width = 1600,
                  height = 700, units = "px")
  
  write.csv(df, paste(results_path, "smajic", analyzed_column, "summary.csv", sep = "_"))
}

## Smajic raw data

for (analyzed_column in column) {
  df <- data.frame(
    n = smajic_raw@meta.data[[analyzed_column]]
  )
  
  table <- summary(smajic_raw@meta.data[[analyzed_column]])
  
  table <- data.frame(
    Statistic = names(table),
    Value = as.numeric(table)
  )
  
  title_path <- paste0("smajic et al., raw cells; n = ", length(df$n))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = n)) + 
    ggplot2::geom_histogram(bins = 40, color = "grey50", fill = "grey78") + 
    ggplot2::labs(title = title_path, 
                  x = analyzed_column, 
                  y = "count") + 
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[4])), color = "black", linetype = "dashed", size = 0.4) + ## Mean
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[3])), color = "#5FCAE3", linetype = "dashed", size = 0.4) +  ## Median
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[2])), color = "#FF9EBC", linetype = "dashed", size = 0.4) + ## 1st Quant
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[5])), color = "#FF9EBC", linetype = "dashed", size = 0.4) +  ## 3rd Quant
    theme_bw()
  
  graph_name <- paste("smajic_raw", analyzed_column,  "summary.png", sep = "_")
  ggplot2::ggsave(graph_name, plot = p,
                  device = "png", path = plot_path, dpi = 200, width = 1600,
                  height = 700, units = "px")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = n)) + 
    ggplot2::geom_histogram(bins = 40, color = "grey50", fill = "grey78") + 
    ggplot2::labs(title = title_path, 
                  x = analyzed_column, 
                  y = "count") + 
    theme_bw()
  
  graph_name <- paste("smajic_raw", analyzed_column,  "summary_without_stats.png", sep = "_")
  ggplot2::ggsave(graph_name, plot = p,
                  device = "png", path = plot_path, dpi = 200, width = 1600,
                  height = 700, units = "px")
  
  write.csv(df, paste(results_path, "smajic_raw", analyzed_column, "summary.csv", sep = "_"))
}

### Kamath

for (analyzed_column in column) {
  df <- data.frame(
    n = kamath@meta.data[[analyzed_column]]
  )
  
  table <- summary(kamath@meta.data[[analyzed_column]])
  
  table <- data.frame(
    Statistic = names(table),
    Value = as.numeric(table)
  )
  
  title_path <- paste0("kamath et al., worked cells; n = ", length(df$n))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = n)) + 
    ggplot2::geom_histogram(bins = 40, color = "grey50", fill = "grey78") + 
    ggplot2::labs(title = title_path, 
                  x = analyzed_column, 
                  y = "count") + 
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[4])), color = "black", linetype = "dashed", size = 0.4) + ## Mean
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[3])), color = "#5FCAE3", linetype = "dashed", size = 0.4) +  ## Median
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[2])), color = "#FF9EBC", linetype = "dashed", size = 0.4) + ## 1st Quant
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[5])), color = "#FF9EBC", linetype = "dashed", size = 0.4) +  ## 3rd Quant
    theme_bw()
  
  graph_name <- paste("kamath", analyzed_column,  "summary.png", sep = "_")
  ggplot2::ggsave(graph_name, plot = p,
                  device = "png", path = plot_path, dpi = 200, width = 1600,
                  height = 700, units = "px")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = n)) + 
    ggplot2::geom_histogram(bins = 40, color = "grey50", fill = "grey78") + 
    ggplot2::labs(title = title_path, 
                  x = analyzed_column, 
                  y = "count") + 
    theme_bw()
  
  graph_name <- paste("kamath", analyzed_column,  "summary_without_stats.png", sep = "_")
  ggplot2::ggsave(graph_name, plot = p,
                  device = "png", path = plot_path, dpi = 200, width = 1600,
                  height = 700, units = "px")
  
  write.csv(df, paste(results_path, "kamath", analyzed_column, "summary.csv", sep = "_"))
}

### Kamath raw 

for (analyzed_column in column) {
  df <- data.frame(
    n = kamath_raw@meta.data[[analyzed_column]]
  )
  
  table <- summary(kamath_raw@meta.data[[analyzed_column]])
  
  table <- data.frame(
    Statistic = names(table),
    Value = as.numeric(table)
  )
  
  title_path <- paste0("kamath et al., raw cells; n = ", length(df$n))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = n)) + 
    ggplot2::geom_histogram(bins = 40, color = "grey50", fill = "grey78") + 
    ggplot2::labs(title = title_path, 
                  x = analyzed_column, 
                  y = "count") + 
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[4])), color = "black", linetype = "dashed", size = 0.4) + ## Mean
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[3])), color = "#5FCAE3", linetype = "dashed", size = 0.4) +  ## Median
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[2])), color = "#FF9EBC", linetype = "dashed", size = 0.4) + ## 1st Quant
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[5])), color = "#FF9EBC", linetype = "dashed", size = 0.4) +  ## 3rd Quant
    theme_bw()
  
  graph_name <- paste("kamath_raw", analyzed_column,  "summary.png", sep = "_")
  ggplot2::ggsave(graph_name, plot = p,
                  device = "png", path = plot_path, dpi = 200, width = 1600,
                  height = 700, units = "px")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = n)) + 
    ggplot2::geom_histogram(bins = 40, color = "grey50", fill = "grey78") + 
    ggplot2::labs(title = title_path, 
                  x = analyzed_column, 
                  y = "count") + 
    theme_bw()
  
  graph_name <- paste("kamath_raw", analyzed_column,  "summary_without_stats.png", sep = "_")
  ggplot2::ggsave(graph_name, plot = p,
                  device = "png", path = plot_path, dpi = 200, width = 1600,
                  height = 700, units = "px")
  
  write.csv(df, paste(results_path, "kamath_raw", analyzed_column, "summary.csv", sep = "_"))
}

### Feleke

for (analyzed_column in column) {
  df <- data.frame(
    n = feleke@meta.data[[analyzed_column]]
  )
  
  table <- summary(feleke@meta.data[[analyzed_column]])
  
  table <- data.frame(
    Statistic = names(table),
    Value = as.numeric(table)
  )
  
  title_path <- paste0("feleke et al., worked cells; n = ", length(df$n))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = n)) + 
    ggplot2::geom_histogram(bins = 40, color = "grey50", fill = "grey78") + 
    ggplot2::labs(title = title_path, 
                  x = analyzed_column, 
                  y = "count") + 
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[4])), color = "black", linetype = "dashed", size = 0.4) + ## Mean
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[3])), color = "#5FCAE3", linetype = "dashed", size = 0.4) +  ## Median
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[2])), color = "#FF9EBC", linetype = "dashed", size = 0.4) + ## 1st Quant
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[5])), color = "#FF9EBC", linetype = "dashed", size = 0.4) +  ## 3rd Quant
    theme_bw()
  
  graph_name <- paste("feleke", analyzed_column,  "summary.png", sep = "_")
  ggplot2::ggsave(graph_name, plot = p,
                  device = "png", path = plot_path, dpi = 200, width = 1600,
                  height = 700, units = "px")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = n)) + 
    ggplot2::geom_histogram(bins = 40, color = "grey50", fill = "grey78") + 
    ggplot2::labs(title = title_path, 
                  x = analyzed_column, 
                  y = "count") + 
    theme_bw()
  
  graph_name <- paste("feleke", analyzed_column,  "summary_without_stats.png", sep = "_")
  ggplot2::ggsave(graph_name, plot = p,
                  device = "png", path = plot_path, dpi = 200, width = 1600,
                  height = 700, units = "px")
  
  write.csv(df, paste(results_path, "feleke", analyzed_column, "summary.csv", sep = "_"))
}

### Feleke raw

for (analyzed_column in column) {
  df <- data.frame(
    n = feleke_raw@meta.data[[analyzed_column]]
  )
  
  table <- summary(feleke_raw@meta.data[[analyzed_column]])
  
  table <- data.frame(
    Statistic = names(table),
    Value = as.numeric(table)
  )
  
  title_path <- paste0("feleke et al., raw cells; n = ", length(df$n))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = n)) + 
    ggplot2::geom_histogram(bins = 40, color = "grey50", fill = "grey78") + 
    ggplot2::labs(title = title_path, 
                  x = analyzed_column, 
                  y = "count") + 
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[4])), color = "black", linetype = "dashed", size = 0.4) + ## Mean
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[3])), color = "#5FCAE3", linetype = "dashed", size = 0.4) +  ## Median
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[2])), color = "#FF9EBC", linetype = "dashed", size = 0.4) + ## 1st Quant
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[5])), color = "#FF9EBC", linetype = "dashed", size = 0.4) +  ## 3rd Quant
    theme_bw()
  
  graph_name <- paste("feleke_raw", analyzed_column,  "summary.png", sep = "_")
  ggplot2::ggsave(graph_name, plot = p,
                  device = "png", path = plot_path, dpi = 200, width = 1600,
                  height = 700, units = "px")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = n)) + 
    ggplot2::geom_histogram(bins = 40, color = "grey50", fill = "grey78") + 
    ggplot2::labs(title = title_path, 
                  x = analyzed_column, 
                  y = "count") + 
    theme_bw()
  
  graph_name <- paste("feleke_raw", analyzed_column,  "summary_without_stats.png", sep = "_")
  ggplot2::ggsave(graph_name, plot = p,
                  device = "png", path = plot_path, dpi = 200, width = 1600,
                  height = 700, units = "px")
  
  write.csv(df, paste(results_path, "feleke_raw", analyzed_column, "summary.csv", sep = "_"))
}

### Yang

for (analyzed_column in column) {
  df <- data.frame(
    n = yang@meta.data[[analyzed_column]]
  )
  
  table <- summary(yang@meta.data[[analyzed_column]])
  
  table <- data.frame(
    Statistic = names(table),
    Value = as.numeric(table)
  )
  
  title_path <- paste0("yang et al., worked cells; n = ", length(df$n))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = n)) + 
    ggplot2::geom_histogram(bins = 40, color = "grey50", fill = "grey78") + 
    ggplot2::labs(title = title_path, 
                  x = analyzed_column, 
                  y = "count") + 
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[4])), color = "black", linetype = "dashed", size = 0.4) + ## Mean
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[3])), color = "#5FCAE3", linetype = "dashed", size = 0.4) +  ## Median
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[2])), color = "#FF9EBC", linetype = "dashed", size = 0.4) + ## 1st Quant
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[5])), color = "#FF9EBC", linetype = "dashed", size = 0.4) +  ## 3rd Quant
    theme_bw()
  
  graph_name <- paste("yang", analyzed_column,  "summary.png", sep = "_")
  ggplot2::ggsave(graph_name, plot = p,
                  device = "png", path = plot_path, dpi = 200, width = 1600,
                  height = 700, units = "px")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = n)) + 
    ggplot2::geom_histogram(bins = 40, color = "grey50", fill = "grey78") + 
    ggplot2::labs(title = title_path, 
                  x = analyzed_column, 
                  y = "count") + 
    theme_bw()
  
  graph_name <- paste("yang", analyzed_column,  "summary_without_stats.png", sep = "_")
  ggplot2::ggsave(graph_name, plot = p,
                  device = "png", path = plot_path, dpi = 200, width = 1600,
                  height = 700, units = "px")
  
  write.csv(df, paste(results_path, "yang", analyzed_column, "summary.csv", sep = "_"))
}

### Yang raw

for (analyzed_column in column) {
  df <- data.frame(
    n = yang_raw@meta.data[[analyzed_column]]
  )
  
  table <- summary(yang_raw@meta.data[[analyzed_column]])
  
  table <- data.frame(
    Statistic = names(table),
    Value = as.numeric(table)
  )
  
  title_path <- paste0("yang et al., raw cells; n = ", length(df$n))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = n)) + 
    ggplot2::geom_histogram(bins = 40, color = "grey50", fill = "grey78") + 
    ggplot2::labs(title = title_path, 
                  x = analyzed_column, 
                  y = "count") + 
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[4])), color = "black", linetype = "dashed", size = 0.4) + ## Mean
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[3])), color = "#5FCAE3", linetype = "dashed", size = 0.4) +  ## Median
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[2])), color = "#FF9EBC", linetype = "dashed", size = 0.4) + ## 1st Quant
    ggplot2::geom_vline(ggplot2::aes(xintercept=as.numeric(table$Value[5])), color = "#FF9EBC", linetype = "dashed", size = 0.4) +  ## 3rd Quant
    theme_bw()
  
  graph_name <- paste("yang_raw", analyzed_column,  "summary.png", sep = "_")
  ggplot2::ggsave(graph_name, plot = p,
                  device = "png", path = plot_path, dpi = 200, width = 1600,
                  height = 700, units = "px")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = n)) + 
    ggplot2::geom_histogram(bins = 40, color = "grey50", fill = "grey78") + 
    ggplot2::labs(title = title_path, 
                  x = analyzed_column, 
                  y = "count") + 
    theme_bw()
  
  graph_name <- paste("yang_raw", analyzed_column,  "summary_without_stats.png", sep = "_")
  ggplot2::ggsave(graph_name, plot = p,
                  device = "png", path = plot_path, dpi = 200, width = 1600,
                  height = 700, units = "px")
  
  write.csv(df, paste(results_path, "yang_raw", analyzed_column, "summary.csv", sep = "_"))
}

### Cell assignation ####

kamath <- readRDS(paste0(our_data, "kamath_preprocessed.rds")) ## Careful, this one is heavy (+/- 36 Gb)
smajic <- readRDS(paste0(our_data, "smajic_preprocessed.rds"))
feleke <- readRDS(paste0(our_data, "feleke_preprocessed.rds"))

assign <- list(Astrocytes = c(3, 16, 26, 34),
               Microglia = c(2, 25),
               Pericytes = c(8),
               SMC = NULL,
               Oligodendrocytes = c(0, 5, 9, 10, 18, 35),
               Endothelial = c(28),
               Neuron = c(1,6,7,12,13,14,15,19,21,22,23,32,33))

kamath@meta.data$named_clusters = "Unknown"
for (clust in names(assign)) {
  if (is.null(assign[[clust]])) {
    next
  } else {
    kamath@meta.data[kamath@meta.data$seurat_clusters %in% assign[[clust]], "named_clusters"] <- clust
  }
}

# Smajic

assign <- list(Astrocytes = c(3, 11, 14), 
               Microglia = c(4), 
               Pericytes = c(9), 
               SMC = NULL, 
               Oligodendrocytes = c(0, 1, 2, 10, 13, 16), 
               Endothelial = c(7), 
               Neuron = c(6,8,12,15,17,18,19))
smajic@meta.data$named_clusters = "Unknown"
for (clust in names(assign)) {
  if (is.null(assign[[clust]])) {
    next
  } else {
    smajic@meta.data[smajic@meta.data$seurat_clusters %in% assign[[clust]], "named_clusters"] <- clust
  }
}

# Feleke

assign <- list(Astrocytes = c(3,5), 
               Microglia = c(13), 
               Pericytes = NULL, 
               SMC = NULL, 
               Oligodendrocytes = c(2, 3, 12), 
               Endothelial = c(19), 
               Neuron = c(0,1,7,9,10,16,17))

feleke@meta.data$named_clusters = "Unknown"
for (clust in names(assign)) {
  if (is.null(assign[[clust]])) {
    next
  } else {
    feleke@meta.data[feleke@meta.data$seurat_clusters %in% assign[[clust]], "named_clusters"] <- clust
  }
}

# Yang 

assign_yang_clusters <- list(Mesenchymal = c("23"), 
                             Endothelial = c("15"), 
                             Ependymal = c("16", "18", "21"), 
                             Epithelial = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                                            "11", "12", "13", "14", "17", "18", "20", "22", "24"), 
                             Macrophage = c("19"))

## Assign a new meta data column to the cell assignation based on the list above. 
yang <- assignation_cluster(yang, assignation_list = assign_yang_clusters, named_column = "approx_clusters", 
                            extraction_cluster = "RNA_snn_res.2")

### Strangely, a base resolution of 0.2 for all three datasets resulted in some datasets not having clusters assign to specific cell types (pericytes and endothelial in kamath and feleke)
## Cluster 28 from Kamath should contain Endothelial
## Cluster 19 from Feleke should contain Pericytes
## Cluster 15 from Yang should contain Mesenchymal cells

### Reclusterisation of Kamath, Feleke, and Yang #### 

kamath_cluster_8 <- subset(kamath, subset = seurat_clusters == "8")
feleke_cluster_19 <- subset(feleke, subset = seurat_clusters == "19")
yang_cluster_15 <- subset(yang, subset = RNA_snn_res.2 == "15")

# Kamath 

setwd(our_data)
SingleCell::analyze_integrated(kamath_cluster_8, organism = "human", assay = "RNA", force_report = TRUE, file_name = "kamath_cluster_8", force_DE = FALSE) 

# Feleke

SingleCell::analyze_integrated(feleke_cluster_19, organism = "human", assay = "RNA", force_report = TRUE, file_name = "feleke_cluster_19")

# Yang

SingleCell::analyze_integrated(yang_cluster_15, organism = "human", assay = "RNA", force_report = TRUE, file_name = "yang_plexus_choroid_control_samples_mensenchymal", force_DE = FALSE,
                               save_path = our_data)

### Cell asssignation 

setwd(path)

kamath_peri <- readRDS(paste0(our_data, "kamath_cluster_8.rds"))
feleke_peri <- readRDS(paste0(our_data, "feleke_cluster_19.rds"))
yang_mesen <- readRDS(paste0(our_data, "yang_plexus_choroid_control_samples_mensenchymal.rds"))


gene_list <- c("GRM8", "PDE7B", "DLC1", "SLC6A12", "TRPC4", "CLDN5", "PECAM1", "VWF")

## This function calls for an equivalent to Seurat::FeaturePlot with the UMAP dimension
series_featureplot(feleke_peri, gene_list, title = "feleke", path = plot_path)
series_featureplot(kamath_peri, gene_list, title = "kamath", path = plot_path)

## Different gene lists for Yang
list_markers_vector <- c("VWF", 
                         "ARL15", "MECOM", "HSP90AA1", "ANO2", 
                         "LINC00276", "HTR2C", "TMEM72-AS1", "GRM8", "PCAT1", "CP", 
                         "TTR", "SOD3", "GPX3", "CTSD", "BSG", "LINGO1", "SLC26A3", "LINC00486", "RAPGEF6", "RASGEF1B",
                         "THSD4", "SLC4A4", "PTGDS", "CEMIP", "KCNMA1",
                         "NRXN1", "DPP10", "CTNNA2", "NPAS3", "FMN2", 
                         "LRMDA", "DOCK4", "ARHGAP24", "SLC8A1", "RUNX1")

translated_gene_list <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = list_markers_vector,
  column = "ENSEMBL", 
  keytype = "SYMBOL"
)

series_featureplot(yang_mesen, translated_gene_list, title = "yang_series", 
                   path = plot_path)

assign_kamath_clusters <- list("8a" = c(0,1,13,14,2,4,6,8,9), 
                               "8b" = c(10,11,12,15,16,17,18,3,5,7))

assign_feleke_clusters <- list("19a" = c(0,11,12,8),
                               "19b" = c(1,10,2,3,4,5,6,7,9))

## Assign a new meta data column to the cell assignation based on the two lists above. 
kamath_peri <- assignation_cluster(kamath_peri, assignation_list = assign_kamath_clusters, named_column = "approx_clusters")
feleke_peri <- assignation_cluster(feleke_peri, assignation_list = assign_feleke_clusters, named_column = "approx_clusters")

assign_kamath <- list(Astrocytes = c(3, 16, 26, 34),
                      Microglia = c(2, 25),
                      Pericytes = c("8a"),
                      SMC = NULL,
                      Oligodendrocytes = c(0, 5, 9, 10, 18, 35),
                      Endothelial = c(28, "8b"),
                      Neuron = c(1,6,7,12,13,14,15,19,21,22,23,32,33))

## assignation_names takes the subsetted seurat, and merge the data with the old seurat object to fuse the new informations with the rest of the data. 
kamath <- assignation_names(kamath_peri, kamath, assignation = assign_kamath)
saveRDS(kamath, paste0(data_path, "kamath_reassembled.rds"))

assign_feleke <- list(Astrocytes = c(3,5), 
                      Microglia = c(13), 
                      Pericytes = c("19a"), 
                      SMC = NULL, 
                      Oligodendrocytes = c(2, 3, 12), 
                      Endothelial = c("19b"), 
                      Neuron = c(0,1,7,9,10,16,17))

feleke <- assignation_names(feleke_peri, feleke, assignation = assign_feleke)
saveRDS(feleke, paste0(data_path, "feleke_reassembled.rds"))

assign_yang_clusters <- list(Mesenchymal = c("0", "1", "2", "3", "4"), 
                             Endothelial = c("5"), 
                             Epithelial = c("6"))

yang_cluster_15 <- assignation_cluster(yang_mesen, assignation_list = assign_yang_clusters, named_column = "approximation_clusters", extraction_cluster = "seurat_clusters")

## Barcodes_reassignation is another assignation function, but this one is based on barcodes instead of another column in meta.data
yang <- barcodes_reassignation(yang, yang_mesen, assignation_list = assign_yang_clusters, new_clusters = "new_clusters", extracted_assignation_subsetted_seurat = "approximation_clusters", 
                               extracted_assignation_seurat = "approx_clusters")

saveRDS(yang, paste0(data_path, "new_assignation_plexus_choroid_yang.rds"))

### Subclusterisation of defined endothelial cells #### 

smajic_endo <- subset(smajic, subset = named_clusters == "Endothelial")
kamath_endo <- subset(kamath, subset = new_assignation == "Endothelial")
feleke_endo <- subset(feleke, subset = new_assignation == "Endothelial")
yang_endo <- subset(yang, subset = new_clusters == "Endothelial")

# Smajic

setwd("~/mnt/scratch/collaborations/JACA/Vasculaire/jaca_vincent/second_analysis/endo_clustering/smajic/")
SingleCell::analyze_integrated(smajic_endo, organism = "human", assay = "RNA", force_report = TRUE, file_name = "smajic_endo.rds", force_DE = FALSE) 

# Kamath 

setwd("~/mnt/scratch/collaborations/JACA/Vasculaire/jaca_vincent/second_analysis/endo_clustering/kamath/")
SingleCell::analyze_integrated(kamath_endo, organism = "human", assay = "RNA", force_report = TRUE, file_name = "kamath_endo.rds", force_DE = FALSE) 

# Feleke

setwd("~/mnt/scratch/collaborations/JACA/Vasculaire/jaca_vincent/second_analysis/endo_clustering/feleke/")
SingleCell::analyze_integrated(feleke_endo, organism = "human", assay = "RNA", force_report = TRUE, file_name = "feleke_endo.rds", force_DE = FALSE)

# Yang

SingleCell::analyze_integrated(yang_endo, organism = "human", assay = "RNA", force_report = TRUE, file_name = "yang_plexus_choroid_control_samples_endothelial.rds", force_DE = FALSE,
                               save_path = "~/mnt/scratch/collaborations/JACA/Vasculaire/jaca_vincent/2025_18_03_Yang_2021/2025_04_10_new_assignation_seurat_object/analyze_integrate_plexus_choroid_control_samples_new_assignation_endothelial/")

### Defining vascular, veinous, and capillary clusters ####

list_markers_smajic <-
  list(
    Arterial = c(
      "ENSG00000150630", "ENSG00000162551" ## VEGFC, ALPL,
    ),
    Cap = c(
      "ENSG00000168389" , "ENSG00000103257" ##  MFSD2A, SLC7A5
    ),
    Veinous = c(
      "ENSG00000115594", "ENSG00000185551" ##  IL1R1, NR2F2
    )
  ) 

list_markers_kamath_feleke <-
  list(
    Arterial = c(
      "VEGFC", "ALPL" 
    ),
    Cap = c(
      "MFSD2A", "SLC7A5" 
    ),
    Veinous = c(
      "IL1R1", "NR2F2"
    )
  ) 

## Smajic

# Exploding the seurat_clusters meta.data into multiple micro clusters
smajic_endo_micro <- Seurat::FindSubCluster(
  smajic_endo,
  cluster = levels(smajic_endo@meta.data$seurat_clusters),
  graph.name = "RNA_snn",
  resolution = 2
)

## cell_assignation_automisation is a function wrapping multiple graphical codes, including a dot plot, which was used based on all 
## micro clusters
cell_assignation_automisation(smajic_endo_micro, list_markers = list_markers_smajic, clusters = "sub.cluster", dotplot_scale = "log.exp", plots_path = plot_path, results_path = results_path, report_path = path, project_name = "GWAS_BEC")

assign <- list(
  Venous = c("0_12", "0_17", "0_2",
             "1_12",  "1_2", 
             "2_12", 
             "3_12", "3_15", "3_18",  
             "4_12", "4_2", 
             "5_12", "5_15", "5_18", "5_2", 
             "6_12", "6_15", "6_2", 
             "7_12", "7_15", "7_17", "7_2"), 
  Capillary = c("0_0", "0_1", "0_10", "0_11", "0_14",  "0_15", "0_4", "0_5", "0_6", "0_7", "0_8",
                "1_0", "1_1", "1_10", "1_13", "1_14", "1_15", "1_17", "1_18", "1_4", "1_5", "1_6", "1_7", "1_8", 
                "2_1", "2_10", "2_11", "2_13", "2_14", "2_17", "2_2", "2_6", "2_7", "2_8", "2_9", 
                "3_1", "3_10", "3_11", "3_13", "3_14", "3_16", "3_17", "3_2", "3_4", "3_6", "3_7", "3_8",
                "4_1", "4_10", "4_11", "4_13", "4_14", "4_4", "4_6", "4_7", "4_8", 
                "5_1", "5_10", "5_11", "5_14", "5_16", "5_4", "5_5", "5_7", "5_8", 
                "6_0", "6_1", "6_10", "6_11", "6_13", "6_14", "6_18", "6_4", "6_6", "6_7", "6_8", "6_9", 
                "7_0", "7_1", "7_10", "7_11", "7_13", "7_14", "7_16", "7_6", "7_7", "7_8"), 
  Arterial = c("0_13", "0_16", "0_3",  "0_9", 
               "1_11",  "1_3",  "1_9", 
               "2_0",  "2_15", "2_16", "2_3", "2_4", "2_5", 
               "3_0",  "3_3", "3_5", "3_9",
               "4_0", "4_15", "4_16", "4_17", "4_3", "4_5", "4_9", 
               "5_0", "5_13", "5_17", "5_3", "5_6", "5_9", 
               "6_16", "6_17", "6_3",  "6_5", 
               "7_18", "7_3", "7_4", "7_5", "7_9")
)

smajic_endo_micro@meta.data$veinous_corrected = "Unknown"
for (clust in names(assign)) {
  if (is.null(assign[[clust]])) {
    next
  } else {
    smajic_endo_micro@meta.data[smajic_endo_micro@meta.data$sub.cluster %in% assign[[clust]], "veinous_corrected"] <- clust
  }
}

## Kamath 

kamath_endo_micro <- Seurat::FindSubCluster(
  kamath_endo,
  cluster = levels(kamath_endo@meta.data$seurat_clusters),
  graph.name = "RNA_snn",
  resolution = 0.5
)

cell_assignation_automisation(kamath_endo_micro, list_markers = list_markers_kamath_feleke, clusters = "sub.cluster", dotplot_scale = "log.exp", plots_path = plot_path, results_path = results_path, report_path = path, project_name = "GWAS_BEC_kamath")

assign <- list(
  Venous = c("0_13", "0_5", "0_6", "0_8", 
             "1_13", "1_3", "1_5", "1_6", "1_8", 
             "10_10", "10_12", "10_5", "10_6", "10_8", 
             "11_13", "11_3", "11_5", "11_6", "11_8", 
             "12_13", "12_5", "12_6", "12_7", "12_8", 
             "13_10", "13_13", "13_5", "13_6", "13_8", 
             "14_13", "14_5", "14_6", "14_8", 
             "15_5", "15_6", "15_8", 
             "16_5", "16_6", "16_8", 
             "17_10", "17_13", "17_5", "17_6", "17_8", 
             "18_5", "18_6", "18_8", 
             "19_10", "19_12", "19_13", "19_5", "19_6", "19_8", 
             "2_13", "2_5", "2_6", "2_8", 
             "3_5", "3_6", "3_8", 
             "4_10", "4_12", "4_13", "4_5", "4_6", "4_8", 
             "5_5", "5_6", "5_8", 
             "6_5", "6_6", "6_8", 
             "7_13", "7_5", "7_6", "7_8",
             "8_5", "8_6", "8_8", 
             "9_10", "9_5", "9_6", "9_8"), 
  
  Capillary = c("0_0", "0_1", "0_10", "0_12", "0_2", "0_3", "0_4", "0_7", "0_9", 
                "1_0", "1_1", "1_10", "1_12", "1_2", "1_4", "1_9", 
                "10_0", "10_1", "10_2", "10_3", "10_4", "10_7", "10_9", 
                "11_0", "11_1", "11_10", "11_12", "11_2", "11_4", "11_7", "11_9", 
                "12_0", "12_1", "12_10", "12_12", "12_2", "12_3", "12_4", "12_9", 
                "13_0", "13_1", "13_12", "13_2", "13_3", "13_4", "13_7", "13_9", 
                "14_0", "14_1", "14_10", "14_12", "14_2", "14_3", "14_4", "14_7", "14_9", 
                "15_0", "15_1", "15_10", "15_12", "15_13", "15_2", "15_3", "15_9", 
                "16_0", "16_1", "16_10", "16_12", "16_2", "16_3", "16_4", "16_7", "16_9", 
                "17_0", "17_1", "17_12", "17_2", "17_3", "17_4", "17_7", "17_9", 
                "18_0", "18_1", "18_10", "18_12", "18_2", "18_3", "18_4", "18_7", "18_9", 
                "19_0", "19_1", "19_2", "19_3", "19_4", "19_7", "19_9", 
                "2_0", "2_1", "2_10", "2_12", "2_2", "2_3", "2_4", "2_7", "2_9", 
                "3_0", "3_1", "3_10", "3_12", "3_13", "3_2", "3_3", "3_4", "3_7", "3_9", 
                "4_0", "4_1", "4_2", "4_3", "4_4", "4_9", 
                "5_0", "5_1", "5_10", "5_12", "5_13", "5_2", "5_3", "5_4", "5_9", 
                "6_0", "6_1", "6_12", "6_13", "6_2", "6_3", "6_4", "6_7", "6_9", 
                "7_0", "7_1", "7_12", "7_2", "7_2", "7_4", "7_7", "7_9", 
                "8_0", "8_1", "8_10", "8_12", "8_13", "8_2", "8_3", "8_4", "8_7", "8_9", 
                "9_0", "9_1", "19_12", "9_13", "9_2", "9_3", "9_4", "9_7", "9_9"), 
  
  Arterial = c("0_11", "0_14", 
               "1_11", "1_7", 
               "10_11", 
               "11_11", "11_14", 
               "12_11", 
               "13_11", 
               "14_11", 
               "15_11", "15_4", "15_7", 
               "16_11", 
               "17_11", 
               "18_11", 
               "19_11", 
               "2_11", 
               "3_11", 
               "4_11", "4_7", 
               "5_11", "5_7", 
               "6_10", "6_11", 
               "7_10", "7_11",
               "8_11", 
               "9_11"),
  
  Unknown = c("10_13")
)

kamath_endo_micro@meta.data$veinous_refined = "Unknown"
for (clust in names(assign)) {
  if (is.null(assign[[clust]])) {
    next
  } else {
    kamath_endo_micro@meta.data[kamath_endo_micro@meta.data$sub.cluster %in% assign[[clust]], "veinous_refined"] <- clust
  }
}

## Feleke

feleke_endo_micro <- Seurat::FindSubCluster(
  feleke_endo,
  cluster = levels(feleke_endo@meta.data$seurat_clusters),
  graph.name = "RNA_snn",
  resolution = 2
)

cell_assignation_automisation(feleke_endo_micro, list_markers = list_markers_kamath_feleke, clusters = "sub.cluster", dotplot_scale = "log.exp", plots_path = plot_path, results_path = results_path, report_path = path, project_name = "GWAS_BEC_feleke")

assign <- list(
  Venous = c("0_0", "0_10", "0_2", "0_9", 
             "1_0", "1_10", "1_2", "1_5", "1_7", "1_9", 
             "2_0", "2_2", "2_7", 
             "3_0", "3_2", "3_9", 
             "4_0", "4_2", "4_4", "4_5", "4_7", 
             "5_0", "5_2", "5_5", "5_7", 
             "6_0", "6_4", "6_7", "6_9", 
             "7_0", "7_1", "7_7", 
             "8_0", "8_4", "8_7", 
             "9_0", "9_2", "9_4"), 
  
  Capillary = c("0_1", "0_3", "0_4", "0_5", 
                "1_1", "1_3", "1_4", "1_6", 
                "2_1", "2_10", "2_3", "2_6", "2_9", 
                "3_1", "3_3", "3_6", "3_7", 
                "4_3", "4_6", "4_9", 
                "5_1", "5_3", "5_4", "5_6", "5_9", 
                "6_1", "6_10", "6_2", "6_3", "6_5", "6_6", 
                "7_10", "7_2", "7_3", "7_4", "7_6", 
                "8_1", "8_2", "8_3", "8_9", 
                "9_1", "9_10", "9_3", "9_5", "9_6", "9_7"), 
  
  Arterial = c("0_8", "1_8", "2_8", "3_4", "3_5", "3_8", 
               "4_1", "4_10", "4_8", "5_8", "6_8", "7_8", 
               "8_6", "8_8", "9_8"), 
  
  Unknown = c("0_6", "5_10", "7_5", "7_9", "8_5", "9_9")
)

feleke_endo_micro@meta.data$veinous_refined = "Unknown"
for (clust in names(assign)) {
  if (is.null(assign[[clust]])) {
    next
  } else {
    feleke_endo_micro@meta.data[feleke_endo_micro@meta.data$sub.cluster %in% assign[[clust]], "veinous_refined"] <- clust
  }
}


### Graphics for Figure 1S -- snRNA-seq datasets captured cell types of the neuro-glia-vascular unit ####

## UMAPs -- Figure 1S.A

pal4 <- c("#1B1A2F", "#4B3B4C", "#E3D581", "#C94C3F", "#EAAE7A" , "#70655A")

# label_clusters_umap is an in-house function that acts like Seurat::FeaturePlot() with a categorical meta.data column. 
label_clusters_umap(smajic, assays = "RNA", slot = "data", path = plot_path, project_name = "smajic", color_palette_ordered = pal4, alpha = 0.6, clusters = "named_clusters", no_borders = TRUE, 
                    width = 900, height = 800)

## Barplots -- Figures 1S.B, 1S.C, and 1S.D

## smajic -- Figures 1S.B ###

cell_counts <- smajic@meta.data %>%
  dplyr::group_by(patient, disease, named_clusters) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::group_by(patient, disease) %>%
  dplyr::mutate(proportion = (n / sum(n))*100) %>%
  dplyr::ungroup() %>%
  dplyr::filter(named_clusters != "Unknown")

stats <- cell_counts %>%
  dplyr::filter(named_clusters != "Unknown") %>%
  dplyr::group_by(named_clusters) %>%
  dplyr::summarise(
    stat = list(wilcox.test(proportion ~ disease, data = cur_data()))
  ) %>%
  dplyr::mutate(
    p_value = sapply(stat, function(x) x$p.value)
  ) %>%
  dplyr::mutate(
    stars = case_when(p_value > 0.05 ~ "NS", 
                      p_value < 0.001 ~ "***", 
                      p_value < 0.01 ~ "**", 
                      p_value < 0.05 ~ "*")
  ) %>%
  dplyr::select(-stat)

y_positions <- cell_counts %>%
  dplyr::group_by(named_clusters) %>%
  dplyr::summarise(y_max = max(proportion)) %>%
  dplyr::left_join(stats, by = "named_clusters") %>%
  dplyr::mutate(y_pos = y_max + 0.05) %>%
  dplyr::filter(stars != "NS")


p <- ggplot2::ggplot(cell_counts, ggplot2::aes(x = named_clusters, y = proportion, fill = disease)) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = position_dodge(width = 0.8),
    color = "black",
    width = 0.7
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.2
  ) + 
  ggplot2::geom_jitter(
    ggplot2::aes(color = disease),
    position = position_jitterdodge(
      jitter.width = 0.15,  
      dodge.width = 0.8     
    ),
    size = 1,
    alpha = 0.5,
    inherit.aes = TRUE
  ) + 
  ggplot2::theme_bw() + 
  ggplot2::labs(
    x = " ",
    y = "% of nuclei",
    fill = "Disease"
  ) + 
  ggplot2::geom_text(data = y_positions,
                     ggplot2::aes(x = named_clusters, y = y_pos, label = stars),
                     inherit.aes = FALSE,
                     vjust = 0, size = 3) +
  ggplot2::scale_fill_manual(values = c("H" = "black", "PD" = "lightsteelblue1")) +
  ggplot2::scale_color_manual(values = c("H" = "grey70", "PD" = "navy")) + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::theme(legend.position = "none")

ggsave(filename = "barplots_proportion_smajic.png", plot = p, path = plot_path,
       dpi = 200, width = 1000, height = 600, units = "px")
write.csv(stats, paste0(results_path, "barplots_statistics_smajic.csv"))
write.csv(cell_counts, paste0(results_path, "barplots_cell_proportion_smajic.csv"))

## kamath -- Figure 1S.C ###

cell_counts <- kamath@meta.data %>%
  dplyr::group_by(patient, disease, new_clusters) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::group_by(patient, disease) %>%
  dplyr::mutate(proportion = (n / sum(n))*100) %>%
  dplyr::ungroup() %>%
  dplyr::filter(new_clusters != "Unknown")

stats <- cell_counts %>%
  dplyr::filter(new_clusters != "Unknown") %>%
  dplyr::group_by(new_clusters) %>%
  dplyr::summarise(
    stat = list(wilcox.test(proportion ~ disease, data = cur_data()))
  ) %>%
  dplyr::mutate(
    p_value = sapply(stat, function(x) x$p.value)
  ) %>%
  dplyr::mutate(
    stars = case_when(p_value > 0.05 ~ "NS", 
                      p_value < 0.001 ~ "***", 
                      p_value < 0.01 ~ "**", 
                      p_value < 0.05 ~ "*")
  ) %>%
  dplyr::select(-stat)

y_positions <- cell_counts %>%
  dplyr::group_by(new_clusters) %>%
  dplyr::summarise(y_max = max(proportion)) %>%
  dplyr::left_join(stats, by = "new_clusters") %>%
  dplyr::mutate(y_pos = y_max + 0.05) %>%
  dplyr::filter(stars != "NS")


p <- ggplot2::ggplot(cell_counts, ggplot2::aes(x = new_clusters, y = proportion, fill = disease)) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = position_dodge(width = 0.8),
    color = "black",
    width = 0.7
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.2
  ) + 
  ggplot2::geom_jitter(
    ggplot2::aes(color = disease),
    position = position_jitterdodge(
      jitter.width = 0.15,  
      dodge.width = 0.8     
    ),
    size = 1,
    alpha = 0.5,
    inherit.aes = TRUE
  ) + 
  ggplot2::theme_bw() + 
  ggplot2::labs(
    x = " ",
    y = "% of nuclei",
    fill = "Disease"
  ) + 
  ggplot2::geom_text(data = y_positions,
                     ggplot2::aes(x = new_clusters, y = y_pos, label = stars),
                     inherit.aes = FALSE,
                     vjust = 0, size = 3) +
  ggplot2::scale_fill_manual(values = c("H" = "black", "PD" = "lightsteelblue1")) +
  ggplot2::scale_color_manual(values = c("H" = "grey70", "PD" = "navy")) + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::theme(legend.position = "none")

p

ggsave(filename = "barplots_proportion_kamath.png", plot = p, path = plot_path,
       dpi = 200, width = 1000, height = 600, units = "px")
write.csv(stats, paste0(results_path, "bbarplots_statistics_kamath.csv"))
write.csv(cell_counts, paste0(results_path, "barplots_cell_proportion_kamath.csv"))

## feleke -- Figure 1S.D ###

cell_counts <- feleke@meta.data %>%
  dplyr::group_by(orig.ident, disease, new_clusters) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::group_by(orig.ident, disease) %>%
  dplyr::mutate(proportion = (n / sum(n))*100) %>%
  dplyr::ungroup() %>%
  dplyr::filter(new_clusters != "Unknown")

stats <- cell_counts %>%
  dplyr::filter(new_clusters != "Unknown") %>%
  dplyr::group_by(new_clusters) %>%
  dplyr::summarise(
    stat = list(wilcox.test(proportion ~ disease, data = cur_data()))
  ) %>%
  dplyr::mutate(
    p_value = sapply(stat, function(x) x$p.value)
  ) %>%
  dplyr::mutate(
    stars = case_when(p_value > 0.05 ~ "NS", 
                      p_value < 0.001 ~ "***", 
                      p_value < 0.01 ~ "**", 
                      p_value < 0.05 ~ "*")
  ) %>%
  dplyr::select(-stat)

y_positions <- cell_counts %>%
  dplyr::group_by(new_clusters) %>%
  dplyr::summarise(y_max = max(proportion)) %>%
  dplyr::left_join(stats, by = "new_clusters") %>%
  dplyr::mutate(y_pos = y_max + 0.05) %>%
  dplyr::filter(stars != "NS")


p <- ggplot2::ggplot(cell_counts, ggplot2::aes(x = new_clusters, y = proportion, fill = disease)) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = position_dodge(width = 0.8),
    color = "black",
    width = 0.7
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.2
  ) + 
  ggplot2::geom_jitter(
    ggplot2::aes(color = disease),
    position = position_jitterdodge(
      jitter.width = 0.15,  
      dodge.width = 0.8     
    ),
    size = 1,
    alpha = 0.5,
    inherit.aes = TRUE
  ) + 
  ggplot2::theme_bw() + 
  ggplot2::labs(
    x = " ",
    y = "% of nuclei",
    fill = "Disease"
  ) + 
  ggplot2::geom_text(data = y_positions,
                     ggplot2::aes(x = new_clusters, y = y_pos, label = stars),
                     inherit.aes = FALSE,
                     vjust = 0, size = 3) +
  ggplot2::scale_fill_manual(values = c("H" = "black", "PD" = "lightsteelblue1")) +
  ggplot2::scale_color_manual(values = c("H" = "grey70", "PD" = "navy")) + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::theme(legend.position = "none")

p

ggsave(filename = "barplots_proportion_feleke.png", plot = p, path = plot_path,
       dpi = 200, width = 1000, height = 600, units = "px")
write.csv(stats, paste0(results_path, "barplots_statistics_feleke.csv"))
write.csv(cell_counts, paste0(results_path, "barplots_cell_proportion_feleke.csv"))

## Barplots by ECs subtypes -- Figures 1S.F, 1S.G, and 1S.H ###

# Smajic Figure 1S.F ###

n_total <- smajic@meta.data %>%
  dplyr::group_by(patient, disease) %>%
  dplyr::summarise(n = n()) 

BEC_counts <- smajic_endo_micro@meta.data %>%
  dplyr::group_by(patient, disease, veinous_corrected) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::left_join(n_total, by = "patient") %>%
  dplyr::mutate(real_proportion = (n.x/n.y) * 100) 

colnames(BEC_counts) <- c("sample", "disease", "veinous", "endothelial_count", "disease_2", "total_cell", "proportion")

table_worked <- BEC_counts %>%
  dplyr::group_by(disease, veinous) %>%
  dplyr::summarise(mean = mean(proportion))

write.csv(BEC_counts, paste0(results_path, "barplots_cell_proportion_smajic_endothelial_proportion.csv"))

p <- ggplot(table_worked, aes(x = disease, y = mean, fill = veinous)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("steelblue", "indianred1", "palegreen2")) +
  labs(x = " ", y = "% of nuclei", fill = "Vascular subtype") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "barplots_proportion_smajic_endothelial.png", plot = p, path = plot_path,
       dpi = 200, width = 900, height = 600, units = "px")

# Kamath Figure 1S.G ###

n_total <- kamath@meta.data %>%
  dplyr::group_by(patient, disease) %>%
  dplyr::summarise(n = n()) 

BEC_counts <- kamath_endo_micro@meta.data %>%
  dplyr::filter(veinous_refined != "Unknown") %>%
  dplyr::group_by(patient, disease, veinous_refined) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::left_join(n_total, by = "patient") %>%
  dplyr::mutate(real_proportion = (n.x/n.y) * 100) 

colnames(BEC_counts) <- c("sample", "disease", "veinous", "endothelial_count", "disease_2", "total_cell", "proportion")

table_worked <- BEC_counts %>%
  dplyr::group_by(disease, veinous) %>%
  dplyr::summarise(mean = mean(proportion))

write.csv(BEC_counts, paste0(results_path, "barplots_cell_proportion_kamath_endothelial_proportion.csv"))

p <- ggplot(table_worked, aes(x = disease, y = mean, fill = veinous)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("steelblue", "indianred1", "palegreen2")) +
  labs(x = " ", y = "% of nuclei", fill = "Vascular subtype") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "barplots_proportion_kamath_endothelial.png", plot = p, path = plot_path,
       dpi = 200, width = 900, height = 600, units = "px")

# Feleke Figure 1S.H ###

n_total <- feleke@meta.data %>%
  dplyr::group_by(orig.ident, disease) %>%
  dplyr::summarise(n = n()) 

BEC_counts <- feleke_endo_micro@meta.data %>%
  dplyr::filter(veinous_refined != "Unknown") %>%
  dplyr::group_by(orig.ident, disease, veinous_refined) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::left_join(n_total, by = "orig.ident") %>%
  dplyr::mutate(real_proportion = (n.x/n.y) * 100) 

colnames(BEC_counts) <- c("sample", "disease", "veinous", "endothelial_count", "disease_2", "total_cell", "proportion")

table_worked <- BEC_counts %>%
  dplyr::group_by(disease, veinous) %>%
  dplyr::summarise(mean = mean(proportion))

write.csv(BEC_counts, paste0(results_path, "barplots_cell_proportion_feleke_endothelial_proportion.csv"))

p <- ggplot(table_worked, aes(x = disease, y = mean, fill = veinous)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("steelblue", "indianred1", "palegreen2")) +
  labs(x = " ", y = "% of nuclei", fill = "Vascular subtype") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "barplots_proportion_feleke_endothelial.png", plot = p, path = plot_path,
       dpi = 200, width = 900, height = 600, units = "px")

### Dotplots -- Graphics for Figure 2S -- Integration of GWAS and snRNA-seq reveals PD risk factors genes enriched at the neuro-glia-vascular unit ####

### Gene lists
gene_endothelial_original <- c("CAVIN2", "ANO2", "ANXA1", "GJA5", "MYLK2", "ELOVL7", "TCIM", "NIBAN2", "CYP1B1", "CD38", "ITIH3", "COL3A1")
gene_microglia_original <- c("CD69", "HLA-DRA", "GPR65", "HLA-DRB5", "FCGR2A", "FCGR2B", "HLA-DQA1", "HLA-DRB1", "STAB1", "MMRN1", "MSR1", "HLA-DQB")
gene_astrocytes_original <- c("CD38", "ITPKB", "APOE", "RFX4", "CASC17", "SLC28A3", "NCRNA00051", "ARD1B", "PBXIP1", "FAM129B", "TRPM2", "LINC00511")
gene_oligodendrocytes_original <- c("DMRT2", "SLC45A3", "PDLIM2", "DNAH17", "LINC00326", "GPM6B", "ARHGEF2", "MAP4K4", "KRTCAP2", "PARP9", "ST6GALNAC3", "C10orf32")

gene_endothelial_feleke <- c("CAVIN2", "ANO2", "ANXA1", "GJA5", "MYLK2", "ELOVL7", "TCIM", "FAM129B", "CYP1B1", "CD38", "ITIH3", "COL3A1") ## FAM129B
gene_microglia_feleke <- c("CD69", "HLA-DRA", "GPR65", "HLA-DRB5", "FCGR2A", "FCGR2B", "HLA-DQA1", "HLA-DRB1", "STAB1", "MMRN1", "MSR1", "HLA-DQB1") ## HLA-DQB1
gene_astrocytes_feleke <- c("CD38", "ITPKB", "APOE", "RFX4", "CASC17", "SLC28A3", "LINC00051", "NAA11", "PBXIP1", "FAM129B", "TRPM2", "LINC00511") ## LINC00051 and NAA11
gene_oligodendrocytes_feleke <- c("DMRT2", "SLC45A3", "PDLIM2", "DNAH17", "LINC00326", "GPM6B", "ARHGEF2", "MAP4K4", "KRTCAP2", "PARP9", "ST6GALNAC3", "BORCS7-ASMT") ## BORCS7-ASMT

gene_endothelial_kamath <- c("SDPR", "ANO2", "ANXA1", "GJA5", "MYLK2", "ELOVL7", "C8orf4", "FAM129B", "CYP1B1", "CD38", "ITIH3", "COL3A1") ## SDPR, FAM129B
gene_microglia_kamath <- c("CD69", "HLA-DRA", "GPR65", "HLA-DRB5", "FCGR2A", "FCGR2B", "HLA-DQA1", "HLA-DRB1", "STAB1", "MMRN1", "MSR1", "HLA-DQB1") ## HLA-DQB1 
gene_astrocytes_kamath  <- c("CD38", "ITPKB", "APOE", "RFX4", "CASC17", "SLC28A3", "LINC00051", "NAA11", "PBXIP1", "FAM129B", "TRPM2", "LINC00511") ## LINC00051 and NAA11
gene_oligodendrocytes_kamath  <- c("DMRT2", "SLC45A3", "PDLIM2", "DNAH17", "LINC00326", "GPM6B", "ARHGEF2", "MAP4K4", "KRTCAP2", "PARP9", "ST6GALNAC3", "C10orf32")

gene_endothelial_smajic <- c("CAVIN2", "ANO2", "ANXA1", "GJA5", "MYLK2", "ELOVL7", "TCIM", "NIBAN2", "CYP1B1", "CD38", "ITIH3", "COL3A1")
gene_microglia_smajic <- c("CD69", "HLA-DRA", "GPR65", "HLA-DRB5", "FCGR2A", "FCGR2B", "HLA-DQA1", "HLA-DRB1", "STAB1", "MMRN1", "MSR1", "HLA-DQB")
gene_astrocytes_smajic <- c("CD38", "ITPKB", "APOE", "RFX4", "CASC17", "SLC28A3", "NCRNA00051", "ARD1B", "PBXIP1", "FAM129B", "TRPM2", "LINC00511")
gene_oligodendrocytes_smajic <- c("DMRT2", "SLC45A3", "PDLIM2", "DNAH17", "LINC00326", "GPM6B", "ARHGEF2", "MAP4K4", "KRTCAP2", "PARP9", "ST6GALNAC3", "C10orf32")

gene_endothelial_translated <- AnnotationDbi::mapIds(
  org.Hs.eg.db, 
  keys = gene_endothelial_smajic,
  column = "ENSEMBL",
  keytype = "ALIAS"
)

gene_microglia_translated <- AnnotationDbi::mapIds(
  org.Hs.eg.db, 
  keys = gene_microglia_smajic,
  column = "ENSEMBL",
  keytype = "ALIAS"
)

gene_astrocytes_translated <- AnnotationDbi::mapIds(
  org.Hs.eg.db, 
  keys = gene_astrocytes_smajic,
  column = "ENSEMBL",
  keytype = "ALIAS"
)

gene_oligodendrocytes_translated <- AnnotationDbi::mapIds(
  org.Hs.eg.db, 
  keys = gene_oligodendrocytes_smajic,
  column = "ENSEMBL",
  keytype = "ALIAS"
)

### The actual dotplots ###

# Smajic ###

## table_dotplot is an inhouse function that recreates Seurat::Dotplot, but only extracts the raw results in a csv table
# This function finds the logaritmn expression of each gene in input (log1p), but also the relative expression, which is used here. 
cell_order <- c("Astrocytes", "Endothelial", "Microglia", "Neuron", "Oligodendrocytes", "Pericytes", "Unknown")
table <- table_dotplot(smajic, gene_order = names(gene_endothelial_translated), genes = gene_endothelial_translated, project_title = "GWAS_smajic_endothelial", id_order = cell_order, taxonomy = "ensembl", clusters = "named_clusters", 
                       path = results_path)

table <- table %>% 
  dplyr::filter(id != "Unknown") ## For each dotplot table, we clean off the "Unknown" cell type.

write.csv(table, paste0("dotplots_tableGWAS_smajic_endothelial.csv"))
## homemade_dotplot is the function that plots the dotplots with the corresponding parameters (either if the info are flipped, if we use the relative expression, color schemes, etc.)
homemade_dotplot(table, gradient_title = "Relative expression", title = "GWAS_smajic_endothelial_list", color = "rel.exp", 
                 path = plot_path, flip = TRUE)

table <- table_dotplot(smajic, gene_order = names(gene_astrocytes_translated), genes = gene_astrocytes_translated, project_title = "GWAS_smajic_astrocytes", id_order = cell_order, taxonomy = "ensembl", clusters = "named_clusters", 
                       path = results_path)

table$symbol[43:49] <- "NCRNA00051" ## Sometimes, the conversion is not perfect with AnnotationDBi, so we have to manually place the 
table$symbol[50:56] <- "ARD1B" ## gene in the table
table$symbol[64:70] <- "FAM129B"
table$symbol[71:77] <- "TRPM2"

table <- table %>% 
  dplyr::filter(id != "Unknown")

write.csv(table, paste0(results_path, "dotplots_tableGWAS_smajic_astrocytes.csv"))
homemade_dotplot(table, gradient_title = "Relative expression", title = "GWAS_smajic_astrocytes_list", color = "rel.exp", 
                 path = plot_path, flip = TRUE)

table <- table_dotplot(smajic, gene_order = names(gene_microglia_translated), genes = gene_microglia_translated, project_title = "GWAS_smajic_microglia", id_order = cell_order, taxonomy = "ensembl", clusters = "named_clusters", 
                       path = results_path)

table$symbol[78:84] <- "HLA-DQB"

table <- table %>% 
  dplyr::filter(id != "Unknown")

write.csv(table, paste0(results_path, "dotplots_tableGWAS_smajic_microglia.csv"))

homemade_dotplot(table, gradient_title = "Relative expression", title = "GWAS_smajic_microglia_list", color = "rel.exp", 
                 path = plot_path, flip = TRUE)

table <- table_dotplot(smajic, gene_order = names(gene_oligodendrocytes_translated), genes = gene_oligodendrocytes_translated, project_title = "GWAS_smajic_oligodendrocytes", id_order = cell_order, taxonomy = "ensembl", clusters = "named_clusters", 
                       path = results_path)

table$symbol[78:84] <- "C10orf32" 
table <- table %>% 
  dplyr::filter(id != "Unknown")

write.csv(table, paste0(results_path, "dotplots_tableGWAS_smajic_oligodendrocytes.csv"))
homemade_dotplot(table, gradient_title = "Relative expression", title = "GWAS_smajic_oligodendrocytes_list", color = "rel.exp", 
                 path = plot_path, flip = TRUE)

## Kamath -- not shown in the final article ###

table <- table_dotplot(kamath, gene_order = gene_endothelial_kamath, genes = gene_endothelial_kamath, project_title = "GWAS_kamath_endothelial", id_order = cell_order, taxonomy = "symbol", clusters = "new_clusters", 
                       path = results_path)

table <- table %>% 
  dplyr::filter(id != "Unknown")

write.csv(table, paste0(results_path, "dotplots_tableGWAS_kamath_endothelial.csv"))
homemade_dotplot(table, gradient_title = "Relative expression", title = "GWAS_kamath_endothelial_list", color = "rel.exp", 
                 path = plot_path, flip = TRUE)

table <- table_dotplot(kamath, gene_order = gene_astrocytes_kamath, genes = gene_astrocytes_kamath, project_title = "GWAS_kamath_astrocytes", id_order = cell_order, taxonomy = "symbol", clusters = "new_clusters", 
                       path = results_path)
table <- table %>% 
  dplyr::filter(id != "Unknown")

write.csv(table, paste0(results_path, "dotplots_tableGWAS_kamath_astrocytes.csv"))
homemade_dotplot(table, gradient_title = "Relative expression", title = "GWAS_kamath_astrocytes_list", color = "rel.exp", 
                 path = plot_path, flip = TRUE)

table <- table_dotplot(kamath, gene_order = gene_microglia_kamath, genes = gene_microglia_kamath, project_title = "GWAS_kamath_microglia", id_order = cell_order, taxonomy = "symbol", clusters = "new_clusters", 
                       path = results_path)

table <- table %>% 
  dplyr::filter(id != "Unknown")

write.csv(table, paste0(results_path, "dotplots_tableGWAS_kamath_microglia.csv"))

homemade_dotplot(table, gradient_title = "Relative expression", title = "GWAS_kamath_microglia_list", color = "rel.exp", 
                 path = plot_path, flip = TRUE)

table <- table_dotplot(kamath, gene_order = gene_oligodendrocytes_kamath, genes = gene_oligodendrocytes_kamath, project_title = "GWAS_kamath_oligodendrocytes", id_order = cell_order, taxonomy = "symbol", clusters = "new_clusters", 
                       path = results_path)

table <- table %>% 
  dplyr::filter(id != "Unknown")

write.csv(table, paste0(results_path, "dotplots_tableGWAS_kamath_oligodendrocytes.csv"))
homemade_dotplot(table, gradient_title = "Relative expression", title = "GWAS_kamath_oligodendrocytes_list", color = "rel.exp", 
                 path = plot_path, flip = TRUE)

## Feleke -- not shown in the article ###

table <- table_dotplot(feleke, gene_order = gene_endothelial_feleke, genes = gene_endothelial_feleke, project_title = "GWAS_feleke_endothelial", id_order = cell_order, taxonomy = "symbol", clusters = "new_clusters", 
                       path = results_path)

table <- table %>% 
  dplyr::filter(id != "Unknown")

write.csv(table, paste0(results_path, "dotplots_tableGWAS_feleke_endothelial.csv"))
homemade_dotplot(table, gradient_title = "Relative expression", title = "GWAS_kamath_endothelial_list", color = "rel.exp", 
                 path = plot_path, flip = TRUE)

table <- table_dotplot(feleke, gene_order = gene_astrocytes_feleke, genes = gene_astrocytes_feleke, project_title = "GWAS_feleke_astrocytes", id_order = cell_order, taxonomy = "symbol", clusters = "new_clusters", 
                       path = results_path)
table <- table %>% 
  dplyr::filter(id != "Unknown")

write.csv(table, paste0(results_path, "dotplots_tableGWAS_feleke_astrocytes.csv"))
homemade_dotplot(table, gradient_title = "Relative expression", title = "GWAS_kamath_astrocytes_list", color = "rel.exp", 
                 path = plot_path, flip = TRUE)

table <- table_dotplot(feleke, gene_order = gene_microglia_feleke, genes = gene_microglia_feleke, project_title = "GWAS_feleke_microglia", id_order = cell_order, taxonomy = "symbol", clusters = "new_clusters", 
                       path = results_path)

table <- table %>% 
  dplyr::filter(id != "Unknown")

write.csv(table, paste0(results_path, "dotplots_tableGWAS_feleke_microglia.csv"))

homemade_dotplot(table, gradient_title = "Relative expression", title = "GWAS_feleke_microglia_list", color = "rel.exp", 
                 path = plot_path, flip = TRUE)

table <- table_dotplot(feleke, gene_order = gene_oligodendrocytes_feleke, genes = gene_oligodendrocytes_feleke, project_title = "GWAS_feleke_oligodendrocytes", id_order = cell_order, taxonomy = "symbol", clusters = "new_clusters", 
                       path = results_path)

table <- table %>% 
  dplyr::filter(id != "Unknown")

write.csv(table, paste0(results_path, "dotplots_tableGWAS_feleke_oligodendrocytes.csv"))
homemade_dotplot(table, gradient_title = "Relative expression", title = "GWAS_feleke_oligodendrocytes_list", color = "rel.exp", 
                 path = plot_path, flip = TRUE)

### Graphics for Figure 3S and Figure 4S (including UMAPs for Figure 1S) ####

gene_list_kamath <- c("SDPR", "ANO2", "ANXA1")
gene_list_feleke <- c("CAVIN2", "ANO2", "ANXA1")
gene_list_smajic <- c("ENSG00000168497", "ENSG00000047617", "ENSG00000135046")
gene_list_smajic_supplement <- c("ENSG00000150630", "ENSG00000162551", "ENSG00000168389", "ENSG00000103257", "ENSG00000115594", "ENSG00000185551")

## Subsetted_UMAPs is an inhouse function that can split a UMAP in two depending on a condition, and append them in one figure

# Figure 1S.E, 3S.D, and 4S.D
subsetted_umaps(smajic_endo, gene_list = gene_list_smajic, taxonomy = "ensembl", title = "GWAS_Smajic", 
                path = plot_path, title_in_italic = TRUE, fixed_limits = TRUE, limit_gradient = c(0,5), no_borders = TRUE)
subsetted_umaps(smajic_endo, gene_list = gene_list_smajic_supplement, taxonomy = "ensembl", title = "GWAS_Smajic", 
                path = plot_path, title_in_italic = TRUE, fixed_limits = TRUE, limit_gradient = c(0,5), no_borders = TRUE)
subsetted_umaps(feleke_endo, gene_list = gene_list_feleke, taxonomy = "taxonomy", title = "GWAS_Feleke", 
                path = plot_path, title_in_italic = TRUE, fixed_limits = TRUE, limit_gradient = c(0,5), no_borders = TRUE)
subsetted_umaps(kamath_endo, gene_list = gene_list_kamath, taxonomy = "taxonomy", title = "GWAS_Kamath", 
                path = plot_path, title_in_italic = TRUE, fixed_limits = TRUE, limit_gradient = c(0,5), no_borders = TRUE)

## Barplots -- Figures 3S.E and 4S.E

# Smajic

n_total <- smajic@meta.data %>%
  dplyr::group_by(patient, disease) %>%
  dplyr::summarise(n = n()) 

BEC_counts <- smajic_endo_micro@meta.data %>%
  dplyr::group_by(patient, disease, veinous_corrected) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::left_join(n_total, by = "patient") %>%
  dplyr::mutate(real_proportion = (n.x/n.y) * 100) 

colnames(BEC_counts) <- c("sample", "disease", "veinous", "endothelial_count", "disease_2", "total_cell", "proportion")

table_worked <- BEC_counts %>%
  dplyr::group_by(disease, veinous) %>%
  dplyr::summarise(mean = mean(proportion))

write.csv(BEC_counts, paste0(results_path, "barplots_cell_proportion_smajic_endothelial_proportion.csv"))

p <- ggplot(table_worked, aes(x = disease, y = mean, fill = veinous)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("steelblue", "indianred1", "palegreen2")) +
  labs(x = " ", y = "% of nuclei", fill = "Vascular subtype") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "barplots_proportion_smajic_endothelial.png", plot = p, path = plot_path,
       dpi = 200, width = 900, height = 600, units = "px")


for (gene in gene_list_smajic) {
  n_total <- smajic@meta.data %>%
    dplyr::group_by(patient, disease) %>%
    dplyr::summarise(n = n()) 
  
  BEC_counts_test <- Seurat::FetchData(smajic_endo_micro, vars = c("patient", "disease", "veinous_corrected", gene)) %>% 
    dplyr::filter(.data[[gene]] > 0.5) %>%
    dplyr::group_by(patient, disease, veinous_corrected) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::left_join(n_total, by = "patient") %>%
    dplyr::mutate(real_proportion = (n.x/n.y) * 100)
  
  colnames(BEC_counts_test) <- c("sample", "disease", "veinous", "endothelial_count", "disease_2", "total_cell", "proportion")
  
  table_worked <- BEC_counts_test %>%
    dplyr::group_by(disease, veinous) %>%
    dplyr::summarise(mean = mean(proportion))
  
  file_name <-  paste0(results_path, "barplots_cell_proportion_smajic_endothelial_", gene, ".csv")
  write.csv(table_worked, file_name)
  
  p <- ggplot2::ggplot(table_worked, ggplot2::aes(x = disease, y = mean, fill = veinous)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::scale_fill_manual(values = c("steelblue", "indianred1", "palegreen2")) +
    ggplot2::labs(x = " ", y = "% of nuclei", fill = "Vascular subtype") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  file_name <- paste0("barplots_proportion_smajic_endothelial_", gene, ".png")
  ggplot2::ggsave(filename = file_name, plot = p, path = plot_path,
                  dpi = 200, width = 900, height = 600, units = "px")
  
}

# Kamath

n_total <- kamath@meta.data %>%
  dplyr::group_by(patient, disease) %>%
  dplyr::summarise(n = n()) 

BEC_counts <- kamath_endo_micro@meta.data %>%
  dplyr::filter(veinous_refined != "Unknown") %>%
  dplyr::group_by(patient, disease, veinous_refined) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::left_join(n_total, by = "patient") %>%
  dplyr::mutate(real_proportion = (n.x/n.y) * 100) 

colnames(BEC_counts) <- c("sample", "disease", "veinous", "endothelial_count", "disease_2", "total_cell", "proportion")

table_worked <- BEC_counts %>%
  dplyr::group_by(disease, veinous) %>%
  dplyr::summarise(mean = mean(proportion))

write.csv(BEC_counts, paste0(results_path, "barplots_cell_proportion_kamath_endothelial_proportion.csv"))

p <- ggplot(table_worked, aes(x = disease, y = mean, fill = veinous)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("steelblue", "indianred1", "palegreen2")) +
  labs(x = " ", y = "% of nuclei", fill = "Vascular subtype") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "barplots_proportion_kamath_endothelial.png", plot = p, path = plot_path,
       dpi = 200, width = 900, height = 600, units = "px")


for (gene in gene_list_kamath) {
  n_total <- kamath@meta.data %>%
    dplyr::group_by(patient, disease) %>%
    dplyr::summarise(n = n()) 
  
  BEC_counts_test <- Seurat::FetchData(kamath_endo_micro, vars = c("patient", "disease", "veinous_refined", gene)) %>% 
    dplyr::filter(.data[[gene]] > 0.5) %>%
    dplyr::filter(veinous_refined != "Unknown") %>%
    dplyr::group_by(patient, disease, veinous_refined) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::left_join(n_total, by = "patient") %>%
    dplyr::mutate(real_proportion = (n.x/n.y) * 100)
  
  colnames(BEC_counts_test) <- c("sample", "disease", "veinous", "endothelial_count", "disease_2", "total_cell", "proportion")
  
  table_worked <- BEC_counts_test %>%
    dplyr::group_by(disease, veinous) %>%
    dplyr::summarise(mean = mean(proportion))
  
  file_name <-  paste0(results_path, "barplots_cell_proportion_kamath_endothelial_", gene, ".csv")
  write.csv(table_worked, file_name)
  
  p <- ggplot2::ggplot(table_worked, ggplot2::aes(x = disease, y = mean, fill = veinous)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::scale_fill_manual(values = c("steelblue", "indianred1", "palegreen2")) +
    ggplot2::labs(x = " ", y = "% of nuclei", fill = "Vascular subtype") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  file_name <- paste0("barplots_proportion_kamath_endothelial_", gene, ".png")
  ggplot2::ggsave(filename = file_name, plot = p, path = plot_path,
                  dpi = 200, width = 900, height = 600, units = "px")
  
}

# Feleke -- not shown in article 

n_total <- feleke@meta.data %>%
  dplyr::group_by(orig.ident, disease) %>%
  dplyr::summarise(n = n()) 

BEC_counts <- feleke_endo_micro@meta.data %>%
  dplyr::filter(veinous_refined != "Unknown") %>%
  dplyr::group_by(orig.ident, disease, veinous_refined) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::left_join(n_total, by = "orig.ident") %>%
  dplyr::mutate(real_proportion = (n.x/n.y) * 100) 

colnames(BEC_counts) <- c("sample", "disease", "veinous", "endothelial_count", "disease_2", "total_cell", "proportion")

table_worked <- BEC_counts %>%
  dplyr::group_by(disease, veinous) %>%
  dplyr::summarise(mean = mean(proportion))

write.csv(BEC_counts, paste0(results_path, "barplots_cell_proportion_feleke_endothelial_proportion.csv"))

p <- ggplot(table_worked, aes(x = disease, y = mean, fill = veinous)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("steelblue", "indianred1", "palegreen2")) +
  labs(x = " ", y = "% of nuclei", fill = "Vascular subtype") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "barplots_proportion_feleke_endothelial.png", plot = p, path = plot_path,
       dpi = 200, width = 900, height = 600, units = "px")


for (gene in gene_list_feleke) {
  n_total <- feleke@meta.data %>%
    dplyr::group_by(orig.ident, disease) %>%
    dplyr::summarise(n = n()) 
  
  BEC_counts_test <- Seurat::FetchData(feleke_endo_micro, vars = c("orig.ident", "disease", "veinous_refined", gene)) %>% 
    dplyr::filter(.data[[gene]] > 0.5) %>%
    dplyr::filter(veinous_refined != "Unknown") %>%
    dplyr::group_by(orig.ident, disease, veinous_refined) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::left_join(n_total, by = "orig.ident") %>%
    dplyr::mutate(real_proportion = (n.x/n.y) * 100)
  
  colnames(BEC_counts_test) <- c("sample", "disease", "veinous", "endothelial_count", "disease_2", "total_cell", "proportion")
  
  table_worked <- BEC_counts_test %>%
    dplyr::group_by(disease, veinous) %>%
    dplyr::summarise(mean = mean(proportion))
  
  file_name <-  paste0(results_path, "barplots_cell_proportion_feleke_endothelial_", gene, ".csv")
  write.csv(table_worked, file_name)
  
  p <- ggplot2::ggplot(table_worked, ggplot2::aes(x = disease, y = mean, fill = veinous)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::scale_fill_manual(values = c("steelblue", "indianred1", "palegreen2")) +
    ggplot2::labs(x = " ", y = "% of nuclei", fill = "Vascular subtype") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  file_name <- paste0("barplots_proportion_feleke_endothelial_", gene, ".png")
  ggplot2::ggsave(filename = file_name, plot = p, path = plot_path,
                  dpi = 200, width = 900, height = 600, units = "px")
  
}

### Graphics for Figure 5S -- Identification of PD risk factors expressed at the BCFSB ####

# Figure 5S.Bi
pal4 <- c("#4B3B4C", "#1B1A2F", "#C94C3F", "#E3D581", "#EAAE7A")

## label_clusters_umap create a UMAP (just like Seurat::FeaturePlot()) but with controlled ggplot2 options -- in house function

label_clusters_umap(yang, assays = "RNA", slot = "data", path = "plots/", project_name = "yang", color_palette_ordered = pal4, alpha = 0.6, clusters = "new_clusters", no_borders = TRUE, 
                    width = 900, height = 800)

# Figure 5S.Bii
yang_count <- Seurat::FetchData(yang, vars = "new_clusters")
yang_count_total <- yang_count %>%
  count(new_clusters)
write_csv(yang_count, file = paste0(results_path, "yang_count.csv"))

x <- ggplot(yang_count_total, aes(x = new_clusters, y = n, fill = new_clusters)) +
  geom_bar(stat = "identity", colour = "grey18") + 
  geom_text(aes(label = n), hjust = -0.120) + 
  theme_bw() + 
  labs(x = "", y = "") + ggplot2::coord_flip() + 
  NoLegend()

ggsave(file = "yang_barplot.png", plot = x, device = "png", 
       path = plot_path,
       dpi = 200, width = 2700, height = 800, units = "px")

### Some analyses were requested for Ependymal cells in Yang ####

yang_epen <- subset(yang, subset = new_clusters == "Ependymal")

data_subsetted <- SingleCell::analyze_integrated(yang_epen, organism = "human", assay = "RNA", force_report = TRUE, file_name = "yang_epen", finding_DEG = FALSE,
                                                 save_path = data_path, step = "normalizing", 
                                                 perform_normalization = FALSE)

gene_list_yang <- c("ENSG00000135046")

df <- data.frame(x = data_subsetted@reductions[["umap"]]@cell.embeddings[,1],
                 y = data_subsetted@reductions[["umap"]]@cell.embeddings[,2],
                 colour_by = Seurat::FetchData(data_subsetted, gene_list_yang))
colnames(df) <- c("x", "y", "colour_by")
max_limit <- ceiling(max(df$colour_by))
min_limit <- floor(min(df$colour_by))
limits <- c(min_limit, max_limit)
first_color = "dodgerblue3"

plot1 <- ggplot2::ggplot(df) + 
  ggplot2::geom_point(ggplot2::aes(x = x, y = y, color = colour_by), alpha = 0.7) + 
  ggplot2::scale_color_gradientn(colours = c("lightgrey", first_color), limits = limits) + 
  ggplot2::labs(title = paste0(gene_list_translated), 
                color = "Expression", 
                x = "umap1", 
                y = "umap2") + 
  ggplot2::theme_bw() + 
  ggplot2::theme(plot.title = element_text(hjust = 0.5))

plot1 <- plot1 + ggplot2::theme_void() + 
  ggplot2::theme(plot.title = element_text(hjust = 0.5))

file_name <- paste0("umaps_ependymal_yang_", gene_list_translated, ".png")
ggplot2::ggsave(file = file_name, plot = plot1, device = "png", 
                path = plot_path,
                dpi = 200, width = 1000, height = 1000, units = "px")


### Statistical analyses (fold change) for genes in Table S1 ####

smajic_no_unknown <- subset(smajic, subset = named_clusters != "Unknown")
kamath_no_unknown <- subset(kamath, subset = new_clusters != "Unknown")
feleke_no_unknown <- subset(feleke, subset = new_clusters != "Unknown")

gene_list_MJFF_GraphA <- c("HSPA1A", "HSP90AA1", "HSPH1", "HSPB1", "DNAJB1", "HSPA1B", "HSPD1", "FKBP4", "HSP90AB1", 
                           "CACYBP", "CHORDC1", "DNAJA1", "BAG3", "PTGES3", "STIP1", "FAM167B", "HSPA8", "CLIC2", 
                           "SNAP23", "IRF1", "AHSA1", "UBB", "ANGPT2", "DNAJB6", "ABHD3", "P4HA1", "IL3RA", "ANKRD37", 
                           "IVNS1ABP", "ST13")

gene_list_MJFF_GraphA_ensembl <- AnnotationDbi::mapIds(
  org.Hs.eg.db, 
  keys = gene_list_MJFF_GraphA, 
  column = "ENSEMBL", 
  keytype = "SYMBOL"
)
gene_list_MJFF_GraphA_ensembl <- unname(gene_list_MJFF_GraphA_ensembl)

# findmarkers_subsetted calculates the fold change on a seurat object. Then a table is simply created as an output.
# This is an equivalent to Seurat::FindMarkers(), but the analyses are made for each category of a meta.data column
# Here, the categorical variable are the cell types.
smajic_table_log <- findmarkers_subsetted(smajic_no_unknown, project_title = "smajic", 
                                          genes = gene_list_MJFF_GraphA_ensembl, path = results_path, taxonomy = "ensembl", clusters = "named_clusters")
kamath_table_log <- findmarkers_subsetted(kamath_no_unknown, project_title = "kamath", 
                                          genes = gene_list_MJFF_GraphA, path = results_path, clusters = "new_clusters")
feleke_table_log <- findmarkers_subsetted(feleke_no_unknown, project_title = "feleke", 
                                          genes = gene_list_MJFF_GraphA, path = results_path, clusters = "new_clusters")

## homemade_dotplot_2FC is an in-house graphical function that takes as input a seurat object, the log_table found previously, and a gene list. 
# This is the same thing as a dot plot. But, because Seurat::Dotplot() cannot create this graph, we had to create a function from scratch. 
# This function is a fancy scatterplot that takes the color gradient as the log2FC at the input table, and the proportion of cells expresseing each gene
# from the seurat object. 
smajic_dotplot <- homemade_dotplot_2FC(smajic, log_table = smajic_table_log, genes = gene_list_MJFF_GraphA_ensembl, project_title = "MJFF_smajic", taxonomy = "ensembl", clusters = "named_clusters",
                                       path = plot_path, flip = TRUE, gene_axis_in_italic = TRUE, 
                                       low_color = "darkblue", mid_color = "seagreen", high_color = "gold", height = 4200) 

kamath_dotplot <- homemade_dotplot_2FC(kamath, log_table = kamath_table_log, clusters = "new_clusters", genes = gene_list_MJFF_GraphA, project_title = "MJFF_kamath", 
                                       path = plot_path, flip = TRUE, gene_axis_in_italic = TRUE, 
                                       low_color = "darkblue", mid_color = "seagreen", high_color = "gold", height = 4200)

feleke_dotplot <- homemade_dotplot_2FC(feleke, log_table = feleke_table_log, clusters = "new_clusters", genes = gene_list_MJFF_GraphA, project_title = "MJFF_feleke", 
                                       path = plot_path, flip = TRUE, gene_axis_in_italic = TRUE, 
                                       low_color = "darkblue", mid_color = "seagreen", high_color = "gold", height = 4200)


### Statistical analyses (fold change) for CAVIN2 and ANXA1 on Kamath, Smajic, and Feleke ####

smajic_list <- c("ENSG00000168497", "ENSG00000135046") 
kamath_list <- c("SDPR", "ANXA1") ## CAVIN2 is not found in Kamath, but this is an equivalent
feleke_list <- c("CAVIN2", "ANXA1")

table <- Seurat::FindMarkers(object = smajic_endo, group.by = "disease", ident.1 = "PD", ident.2 = "H", min.pct = 0.05, logfc.threshold = 0, 
                             test.use = "wilcox", features = smajic_list)

rownames(table) <- c("CAVIN2", "ANXA1")
write.csv(table, paste0(results_path, "DEG_analyses_smajic_endo.csv"))

table <- Seurat::FindMarkers(object = kamath_endo, group.by = "disease", ident.1 = "PD", ident.2 = "H", min.pct = 0.05, logfc.threshold = 0, 
                             test.use = "wilcox", features = kamath_list)

write.csv(table, paste0(results_path, "DEG_analyses_kamath_endo.csv"))


table <- Seurat::FindMarkers(object = feleke_endo, group.by = "disease", ident.1 = "PD", ident.2 = "H", min.pct = 0.05, logfc.threshold = 0, 
                             test.use = "wilcox", features = feleke_list)

write.csv(table, paste0(results_path, "DEG_analyses_feleke_endo.csv")) ## Non significant results

## Make the same analyses on MBECs data ### 

MBECs <- c("Venous", "Capillary", "Arterial")

# Smajic ##

for (bec in MBECs) {
  
  subsetted <- subset(smajic_endo_micro, subset = veinous_corrected == bec)
  
  table <- Seurat::FindMarkers(object = smajic_endo_micro, group.by = "disease", ident.1 = "PD", ident.2 = "H", 
                               min.pct = 0.05, logfc.threshold = 0, 
                               test.use = "wilcox", features = smajic_list)
  
  table_name <- paste0(results_path, "DEG_analyses_smajic_endo_", bec, ".csv")
  write.csv(table, table_name)
}

# Kamath ## 

for (bec in MBECs) {
  
  subsetted <- subset(kamath_endo_micro, subset = veinous_refined == bec)
  
  table <- Seurat::FindMarkers(object = kamath_endo_micro, group.by = "disease", ident.1 = "PD", ident.2 = "H", 
                               min.pct = 0.05, logfc.threshold = 0, 
                               test.use = "wilcox", features = kamath_list)
  
  table_name <- paste0(results_path, "DEG_analyses_kamath_endo_", bec, ".csv")
  write.csv(table, table_name)
}

# Feleke ## 

for (bec in MBECs) {
  
  subsetted <- subset(feleke_endo_micro, subset = veinous_refined == bec)
  
  table <- Seurat::FindMarkers(object = feleke_endo_micro, group.by = "disease", ident.1 = "PD", ident.2 = "H", 
                               min.pct = 0.05, logfc.threshold = 0, 
                               test.use = "wilcox", features = feleke_list)
  
  table_name <- paste0(results_path, "DEG_analyses_kamath_endo_", bec, ".csv")
  write.csv(table, table_name)
}

### Statistical analyses (fold change) for some genes on neurons for Kamath, Smajic, and Feleke ####

gene_list_smajic <- c("ENSG00000180176", "ENSG00000165092", "ENSG00000110693") ## PITX3 is not found
gene_list_kamath_feleke <- c("TH", "PITX3", "ALDH1A1", "SOX6")

kamath_neuro <- subset(kamath, subset = new_clusters == "Neuron")
smajic_neuro <- subset(smajic, subset = named_clusters == "Neuron")
feleke_neuro <- subset(feleke, subset = new_clusters == "Neuron")

table_smajic <- Seurat::FindMarkers(object = smajic_neuro, group.by = "disease", ident.1 = "PD", ident.2 = "H", min.pct = 0, logfc.threshold = 0, 
                             test.use = "wilcox", features = gene_list_smajic)

write.csv(table_smajic, paste0(results_path, "DEG_analyses_smajic_neuro.csv"))

table <- Seurat::FindMarkers(object = kamath_neuro, group.by = "disease", ident.1 = "PD", ident.2 = "H", min.pct = 0, logfc.threshold = 0, 
                             test.use = "wilcox", features = gene_list_kamath_feleke)

write.csv(table, paste0(results_path, "DEG_analyses_kamath_neuro.csv"))


table <- Seurat::FindMarkers(object = feleke_neuro, group.by = "disease", ident.1 = "PD", ident.2 = "H", min.pct = 0, logfc.threshold = 0, 
                             test.use = "wilcox", features = gene_list_kamath_feleke)

write.csv(table, paste0(results_path, "DEG_analyses_feleke_neuro.csv"))

### Making the GWAS list and the heatmaps per datasets #### 

gwas <- read_tsv(paste0(our_data, "gwas.tsv")) ## From the https://www.ebi.ac.uk/gwas/docs/file-downloads catalog from GWAS (v1.0.2.1)

gwas <- gwas %>% dplyr::mutate(across(c(`REPORTED GENE(S)`, MAPPED_GENE), ~ gsub(pattern = " - ", replacement = ", ", x = .))) %>% separate_longer_delim(`REPORTED GENE(S)`, delim = ", ") %>% separate_longer_delim(MAPPED_GENE, delim = ", ")
gwas_genes <- unique(c(gwas$`REPORTED GENE(S)`, gwas$MAPPED_GENE))

# smajic
i = "smajic"
scdata <- readRDS(paste0(our_data, i, "_preprocessed.rds"))

assign <- list(Astrocytes = c(3, 11, 14), 
               Microglia = c(4), 
               Pericytes = c(9), 
               SMC = NULL, 
               Oligodendrocytes = c(0, 1, 2, 10, 13, 16), 
               Endothelial = c(7), 
               Neuron = c(6,8,12,15,17,18,19))

scdata@meta.data$named_clusters = "Unknown"
for (clust in names(assign)) {
  if (is.null(assign[[clust]])) {
    next
  } else {
    scdata@meta.data[scdata@meta.data$seurat_clusters %in% assign[[clust]], "named_clusters"] <- clust
  }
}

if (i %in% c("smajic")) {
  genes <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = gwas_genes,
    column = "ENSEMBL",
    keytype = "ALIAS")
} else {
  genes <- gwas_genes
}
Idents(scdata) <- "named_clusters"

genes2 <- genes[!is.na(genes)]

avg <- as.data.frame(Seurat::AverageExpression(scdata))
not_genes <- c(genes2[!genes2 %in% row.names(avg)], genes[is.na(genes)])
genes <- genes2[genes2 %in% row.names(avg)]
avg <- avg[unlist(genes),]
row.names(avg) <- names(genes)
colnames(avg) <- gsub("RNA.", "", colnames(avg))

write.csv(avg, paste0(results_path, i, "_gwas_mean_expression.csv"))

pdf(file = paste0(plot_path, i, "_gwas_genes_heatmap.pdf"), width = 5, height = 12)
hm <- ComplexHeatmap::Heatmap(avg,
                              column_title = paste0("gwas genes for dataset ", i),
                              show_row_names = FALSE,
                              cluster_columns = FALSE,
                              cluster_rows = TRUE)
hm <- draw(hm)
dev.off()

### feleke
i = "feleke"
scdata <- readRDS(paste0(our_data, i, "_preprocessed.rds"))

assign <- list(Astrocytes = c(3,5), 
               Microglia = c(13), 
               Pericytes = NULL, 
               SMC = NULL, 
               Oligodendrocytes = c(2, 3, 12), 
               Endothelial = c(19), 
               Neuron = c(0,1,7,9,10,16,17))

scdata@meta.data$named_clusters = "Unknown"
for (clust in names(assign)) {
  if (is.null(assign[[clust]])) {
    next
  } else {
    scdata@meta.data[scdata@meta.data$seurat_clusters %in% assign[[clust]], "named_clusters"] <- clust
  }
}

if (i %in% c("smajic")) {
  genes <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = gwas_genes,
    column = "ENSEMBL",
    keytype = "ALIAS")
} else {
  genes <- gwas_genes
}
Idents(scdata) <- "named_clusters"

genes2 <- genes[!is.na(genes)]

avg <- as.data.frame(Seurat::AverageExpression(scdata))
not_genes <- c(genes2[!genes2 %in% row.names(avg)], genes[is.na(genes)])
genes <- genes2[genes2 %in% row.names(avg)]
avg <- avg[unlist(genes),]

colnames(avg) <- gsub("RNA.", "", colnames(avg))
write.csv(avg, paste0(results_path, i, "_gwas_mean_expression.csv"))

pdf(file = paste0(plot_path, i, "_gwas_genes_heatmap.pdf"), width = 5, height = 12)
hm <- ComplexHeatmap::Heatmap(avg,
                              column_title = paste0("gwas genes for dataset ", i),
                              show_row_names = FALSE,
                              cluster_columns = FALSE,
                              cluster_rows = TRUE) 
hm <- draw(hm)
dev.off()


### kamath
i = "kamath"
scdata <- readRDS(paste0(our_data, i, "_preprocessed.rds"))

assign <- list(Astrocytes = c(3, 16, 26, 34),
               Microglia = c(2, 25),
               Pericytes = c(8),
               SMC = NULL,
               Oligodendrocytes = c(0, 5, 9, 10, 18, 35),
               Endothelial = c(28),
               Neuron = c(1,6,7,12,13,14,15,19,21,22,23,32,33))

scdata@meta.data$named_clusters = "Unknown"
for (clust in names(assign)) {
  if (is.null(assign[[clust]])) {
    next
  } else {
    scdata@meta.data[scdata@meta.data$seurat_clusters %in% assign[[clust]], "named_clusters"] <- clust
  }
}
Idents(scdata) <- "named_clusters"

genes <- gwas_genes

genes2 <- genes[!is.na(genes)]

avg <- as.data.frame(Seurat::AverageExpression(scdata))
not_genes <- c(genes2[!genes2 %in% row.names(avg)], genes[is.na(genes)])
genes <- genes2[genes2 %in% row.names(avg)]
avg <- avg[unlist(genes),]

colnames(avg) <- gsub("RNA.", "", colnames(avg))
write.csv(avg, paste0(results_path, i, "_gwas_mean_expression.csv"))

### Average Expression tables in excel ####

# Gene lists 

### Files worked on the previous section
kamath_table_gwas <- read.csv(paste0(results_path, "kamath_gwas_mean_expression.csv"))
smajic_table_gwas <- read.csv(paste0(results_path, "smajic_gwas_mean_expression.csv"))
feleke_table_gwas <- read.csv(paste0(results_path, "feleke_gwas_mean_expression.csv"))

kamath_list <- kamath_table_gwas$X
smajic_list <- smajic_table_gwas$X
feleke_list <- feleke_table_gwas$X

order <- c("Neuron", "Endothelial", "Pericytes", "Microglia", "Astrocytes", "Oligodendrocytes")

## average_expression_annotation wraps the Seurat::AverageExpression() function from Seurat, but can splitted the data based on the meta.data selected, splitted the data based on a treatment, 
# and arrange the tables in function of a defined order (object "order" in this case)

smajic_PD <- average_expression_annotation(smajic, project_title = "smajic", gene_list = smajic_list,  subset_treatment = "PD", taxonomy = "ensembl", order = order, skip_unknown = TRUE)
smajic_H <- average_expression_annotation(smajic, project_title = "smajic", gene_list = smajic_list,  subset_treatment = "H", taxonomy = "ensembl", order = order, skip_unknown = TRUE)
kamath_PD <- average_expression_annotation(kamath, project_title = "kamath", gene_list = kamath_list, group.by = "new_assignation",  subset_treatment = "PD", taxonomy = "symbol", order = order, skip_unknown = TRUE)
kamath_H <- average_expression_annotation(kamath, project_title = "kamath", gene_list = kamath_list,  group.by = "new_assignation", subset_treatment = "H", taxonomy = "symbol", order = order, skip_unknown = TRUE)
feleke_PD <- average_expression_annotation(feleke, project_title = "feleke", gene_list = feleke_list,  group.by = "new_assignation", subset_treatment = "PD", taxonomy = "symbol", order = order, skip_unknown = TRUE)
feleke_H <- average_expression_annotation(feleke, project_title = "feleke", gene_list = feleke_list, group.by = "new_assignation", subset_treatment = "H", taxonomy = "symbol", order = order, skip_unknown = TRUE)

merged_table_PD <- merge(feleke_PD, kamath_PD, by = "gene_name", all = TRUE)
merged_table_PD <- merge(merged_table_PD, smajic_PD, by = "gene_name", all = TRUE)
merged_table_PD <- arrange(merged_table_PD, gene_name)
merged_table_PD <- relocate(merged_table_PD, 
                            gene_name,
                            Neuron_smajic_PD, Neuron_feleke_PD, Neuron_kamath_PD,
                            Endothelial_smajic_PD, Endothelial_feleke_PD, Endothelial_kamath_PD,
                            Pericytes_smajic_PD, Pericytes_feleke_PD, Pericytes_kamath_PD,
                            Astrocytes_smajic_PD, Astrocytes_feleke_PD, Astrocytes_kamath_PD,
                            Microglia_smajic_PD, Microglia_feleke_PD, Microglia_kamath_PD,  
                            Oligodendrocytes_smajic_PD, Oligodendrocytes_feleke_PD, Oligodendrocytes_kamath_PD,  
)

merged_table_H <- merge(feleke_H, kamath_H, by = "gene_name", all = TRUE)
merged_table_H <- merge(merged_table_H, smajic_H, by = "gene_name", all = TRUE)
merged_table_H <- arrange(merged_table_H, gene_name)
merged_table_H <- relocate(merged_table_H, 
                           gene_name,
                           Neuron_smajic_H, Neuron_feleke_H, Neuron_kamath_H,
                           Endothelial_smajic_H, Endothelial_feleke_H, Endothelial_kamath_H,
                           Pericytes_smajic_H, Pericytes_feleke_H, Pericytes_kamath_H,
                           Astrocytes_smajic_H, Astrocytes_feleke_H, Astrocytes_kamath_H,
                           Microglia_smajic_H, Microglia_feleke_H, Microglia_kamath_H,  
                           Oligodendrocytes_smajic_H, Oligodendrocytes_feleke_H, Oligodendrocytes_kamath_H,  
)

# Mean per cell types 

mean_table_PD <- merged_table_PD %>%
  mutate(
    Neuron_PD = rowMeans(cbind(Neuron_smajic_PD, Neuron_feleke_PD, Neuron_kamath_PD), na.rm = TRUE),
    Endothelial_PD = rowMeans(cbind(Endothelial_smajic_PD, Endothelial_feleke_PD, Endothelial_kamath_PD), na.rm = TRUE),
    Pericytes_PD = rowMeans(cbind(Pericytes_smajic_PD, Pericytes_feleke_PD, Pericytes_kamath_PD), na.rm = TRUE),
    Astrocytes_PD = rowMeans(cbind(Astrocytes_smajic_PD, Astrocytes_feleke_PD, Astrocytes_kamath_PD), na.rm = TRUE),
    Microglia_PD = rowMeans(cbind(Microglia_smajic_PD, Microglia_feleke_PD, Microglia_kamath_PD), na.rm = TRUE),
    Oligodendrocytes_PD = rowMeans(cbind(Oligodendrocytes_smajic_PD, Oligodendrocytes_feleke_PD, Oligodendrocytes_kamath_PD), na.rm = TRUE)
  )

mean_table_PD <- mean_table_PD[, c(1, 20, 21, 22, 23, 24, 25)]

mean_table_H <- merged_table_H %>%
  mutate(
    Neuron_H = rowMeans(cbind(Neuron_smajic_H, Neuron_feleke_H, Neuron_kamath_H), na.rm = TRUE),
    Endothelial_H = rowMeans(cbind(Endothelial_smajic_H, Endothelial_feleke_H, Endothelial_kamath_H), na.rm = TRUE),
    Pericytes_H = rowMeans(cbind(Pericytes_smajic_H, Pericytes_feleke_H, Pericytes_kamath_H), na.rm = TRUE),
    Astrocytes_H = rowMeans(cbind(Astrocytes_smajic_H, Astrocytes_feleke_H, Astrocytes_kamath_H), na.rm = TRUE),
    Microglia_H = rowMeans(cbind(Microglia_smajic_H, Microglia_feleke_H, Microglia_kamath_H), na.rm = TRUE),
    Oligodendrocytes_H = rowMeans(cbind(Oligodendrocytes_smajic_H, Oligodendrocytes_feleke_H, Oligodendrocytes_kamath_H), na.rm = TRUE)
  )

mean_table_H <- mean_table_H[, c(1, 20, 21, 22, 23, 24, 25)]

mean_max_table_PD <- mean_table_PD
str(mean_max_table_PD)
mean_max_table_PD$Max_PD <- apply(mean_max_table_PD[, -1], 1, function(x) max(x, na.rm = TRUE)) 
mean_max_table_PD$Max_PD_Cell_Type <- apply(mean_max_table_PD[, -1], 1, function(x) names(x)[which.max(x)])

mean_max_table_H <- mean_table_H
str(mean_max_table_H)
mean_max_table_H$Max_H <- apply(mean_max_table_H[, -1], 1, function(x) max(x, na.rm = TRUE)) 
mean_max_table_H$Max_H_Cell_Type <- apply(mean_max_table_H[, -1], 1, function(x) names(x)[which.max(x)])

table_list <- list(
  PD_merged <- merged_table_PD,
  PD_mean <- mean_max_table_PD,
  H_merged <- merged_table_H,
  H_mean <- mean_max_table_H
)
sheet_names <- c("PD_merged", "PD_mean","H_merged", "H_mean")

openxlsx::write.xlsx(x = table_list, sheetName = sheet_names, 
                     file = paste0(results_path, "GWAS_AverageExpression_new_list.xlsx"))

gene_names_GWAS <- mean_table_H[, 1]
order <- c("Endothelial", "Ependymal", "Epithelial", "Macrophage", "Mesenchymal")

## For Yang et al., specifically 

yang_H <- average_expression_annotation(yang, project_title = "yang", gene_list = gene_names_GWAS,  taxonomy = "ensembl", group.by = "new_clusters", order = order, skip_unknown = FALSE, treatment = "orig.ident") 

table_list <- list(
  PD_merged <- merged_table_PD,
  PD_mean <- mean_max_table_PD,
  H_merged <- merged_table_H,
  H_mean <- mean_max_table_H, 
  Yang_mean <- yang_H
)

sheet_names <- c("PD_merged", "PD_mean","H_merged", "H_mean", "Yang_mean")

openxlsx::write.xlsx(x = table_list, sheetName = sheet_names, 
                     file = paste0(results_path, "GWAS_AverageExpression_new_list.xlsx"))