#!/usr/bin/env Rscript

# START ----
# Workflow for processing combined analysis object
#  - Prior to running this code set working directory to 'impact-analysis/code'

library(Seurat)
library(SeuratData)
library(rjson)
library(tidyverse)
library(stringr)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(SingleR)
library(Azimuth)
library(kableExtra)
library(devtools)
library(parallel)
library(data.table)
library(harmony)
library(DoubletFinder)
library(impactSingleCellToolkit)

# INPUT and OUTPUT Directories
# 2_preprocessing_combined_object_gex.R  [Param File]

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("ERROR: At least one argument must be supplied (JSON parameter file).json", call.=FALSE)
} 
param_file_fh = args[1]

#param_file_fh = "../input/input_one_patient_analysis.json"
params        = fromJSON(file = param_file_fh)

# INPUT
metadata_fh   = params$metadata_file
metadata      = read.csv(metadata_fh,stringsAsFactors = F,check.names = F)
THREADS       = params$threads
PERSAMPLE_SCT = FALSE

# OUTPUT Results Directories
plot_dir     = params$plot_directory
objects_dir  = params$objects_directory
if(!file.exists(plot_dir)){dir.create(plot_dir,recursive = TRUE)}
if(!file.exists(objects_dir)){dir.create(objects_dir,recursive = TRUE)}

# Input and Output RData Objects
fh_raw_seurat_obj       <- file.path(objects_dir,"seurat_obj.combined.gex.RData")
fh_processed_seurat_obj <- file.path(objects_dir,"seurat_obj.processed.rds")

# Increase size of ENV
options(future.globals.maxSize= 891289600)
set.seed(112358)



# 1. Load and start normalization/clustering ----

# Load combined GEX object
load(fh_raw_seurat_obj)

# Normalization per sample
if(PERSAMPLE_SCT){
  seurat.obj.list <- SplitObject(seurat.obj, split.by="sample")
  seurat.obj.list <- lapply(X = seurat.obj.list, 
                            FUN = SCTransform, 
                            return.only.var.genes = FALSE)
  seurat.obj <- merge(seurat.obj.list[[1]], 
                      y = seurat.obj.list[2:length(seurat.obj.list)], 
                      add.cell.ids = names(seurat.obj.list))
  DefaultAssay(seurat.obj) <- "SCT"
}


# Re-cluster
seurat.obj   <- scaleAndClusterSeuratObject(seurat.obj,normalize = T,dims = 1:30,npca = 10,tsne = T)
pc_elbowplot <- plotOptimalPCsforSeuratObject(seurat.obj)

# 2. Plot pre-batch correction and add in batch column ----

# add batch column
metadata_gex                <- metadata[metadata$type %in% c("counts_gex","counts"),]
sample_to_batch_list        <- metadata_gex$run
names(sample_to_batch_list) <- metadata_gex$sample
sample_to_batch_list        <- as.list(sample_to_batch_list)

batch_col   <- seurat.obj$sample
uni_samples <- unique(batch_col)
for (s in uni_samples) {
  batch <- sample_to_batch_list[[s]]
  batch_col[batch_col %in% s] <- batch
}
seurat.obj$batch <- batch_col

# Plot 
plot_pre_batch_correction <- UMAPPlot(object = seurat.obj, group.by="batch")

png(file.path(plot_dir,"pre_batch_correction_umap.png"),height = 500,width = 600)
print(plot_pre_batch_correction)
dev.off()

# 3. Harmony batch correction ----

seurat.obj <- seurat.obj %>% RunHarmony("batch", assay.use="SCT")

# UMAP and clustering with harmonized PCs
seurat.obj <- RunUMAP(seurat.obj, reduction='harmony', dims = 1:30)
seurat.obj <- FindNeighbors(seurat.obj, reduction='harmony')
seurat.obj <- FindClusters(seurat.obj, resolution = 0.3)

# Plot 
plot_post_batch_correction <- UMAPPlot(object = seurat.obj, group.by="batch")

png(file.path(plot_dir,"post_batch_correction_umap.png"),height = 500,width = 600)
print(plot_post_batch_correction)
dev.off()


# 4. Azimuth automated celltype annotations ----
###InstallData("pbmcref")
###DefaultAssay(seurat.obj). # Check that assay is either RNA or SCT, not CSP

# Run Azimuth
results_azimuth  <- RunAzimuth(seurat.obj, reference = "pbmcref")

# Add annotations back into Seurat object
celltypes_azimuth <- results_azimuth@meta.data
annot_cols        <- names(celltypes_azimuth)
annot_cols        <- annot_cols[annot_cols %in% c("predicted.celltype.l1","predicted.celltype.l2","predicted.celltype.l3")]
celltypes_azimuth <- celltypes_azimuth[,annot_cols]

seurat.obj@meta.data <- cbind(seurat.obj@meta.data,celltypes_azimuth)


# 5. UMAP celltype annotations ----

# Plot Azimuth annotations

### Level 1 annotations (main)
l1_cell_annot_azimuth <- UMAPPlot(object = seurat.obj, group.by="predicted.celltype.l1")

png(file.path(plot_dir,"celltype_annot_azimuth_l1_umap.png"),height = 500,width = 600)
print(l1_cell_annot_azimuth)
dev.off()

### Level 2 annotations (fine)
l2_cell_annot_azimuth <- UMAPPlot(object = seurat.obj, group.by="predicted.celltype.l2")

png(file.path(plot_dir,"celltype_annot_azimuth_l2_umap.png"),height = 500,width = 800)
print(l2_cell_annot_azimuth)
dev.off()


# 6. Process CSP part of object ----

# This section is option depending on if Cite-Seq data are added in the object
if("CSP" %in% names(seurat.obj)){
  DefaultAssay(seurat.obj) <- "CSP" # it should be 'SCT' after processing the GEX data
  
  ### Normalize
  seurat.obj <- NormalizeData(seurat.obj, normalization.method = "CLR", margin = 2, assay = "CSP")
  ### Run PCA
  VariableFeatures(seurat.obj) = rownames(seurat.obj@assays[["CSP"]])
  seurat.obj = seurat.obj %>% 
    ScaleData(verbose=F) %>%
    RunPCA(reduction.name="apca",approx=F, verbose=F) 
  ### Pick PCs
  total_variance <- sum(matrixStats::rowVars(
    as.matrix(seurat.obj@assays[["CSP"]]@scale.data)))
  eigValues = (seurat.obj@reductions$apca@stdev)^2  
  varExplained = eigValues / total_variance
  csp_pc_plot <- plot(varExplained)
  
  ### Run UMAP
  seurat.obj <- RunUMAP(seurat.obj, 
                        reduction = 'apca', 
                        dims = 1:12, 
                        assay = 'CSP', 
                        reduction.name = 'csp.umap', 
                        reduction.key = 'cspUMAP_',
                        verbose = F)
}


# SAVE ----
DefaultAssay(seurat.obj) <- "SCT"
saveRDS(seurat.obj,file = fh_processed_seurat_obj)
