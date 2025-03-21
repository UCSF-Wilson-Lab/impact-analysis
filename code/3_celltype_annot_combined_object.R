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
# 3_celltype_annot_combined_object.R  [Param File]

#args = commandArgs(trailingOnly=TRUE)
#if (length(args)!=1) {
#  stop("ERROR: At least one argument must be supplied (JSON parameter file).json", call.=FALSE)
#} 
#param_file_fh = args[1]

#param_file_fh = "../input/input_one_patient_analysis.json"
param_file_fh = "/wynton/protected/home/wilson/rdandekar/rprojects/ImpactAnalysis/impact-analysis/input/input_all_csf_analysis.json"
params        = fromJSON(file = param_file_fh)

# INPUT
metadata_fh       = params$metadata_file
metadata          = read.csv(metadata_fh,stringsAsFactors = F,check.names = F)
ref_dir_azimuth   = params$celltype_annotation_reference_directory
THREADS           = params$threads

# OUTPUT Results Directories
plot_dir     = params$plot_directory
objects_dir  = params$objects_directory
if(!file.exists(plot_dir)){dir.create(plot_dir,recursive = TRUE)}
if(!file.exists(objects_dir)){dir.create(objects_dir,recursive = TRUE)}


# Input and Output RData Objects
fh_raw_seurat_obj       <- file.path(objects_dir,"seurat_obj.combined.gex.TEMP.RData")
fh_processed_seurat_obj <- file.path(objects_dir,"seurat_obj.processed.rds")

# Increase size of ENV
options(future.globals.maxSize= 40000*1024^2)
set.seed(112358)



# 1. Load and start normalization/clustering ----

# Load combined GEX object
load(fh_raw_seurat_obj)


# 2. Azimuth automated celltype annotations ----
###InstallData("pbmcref")
###DefaultAssay(seurat.obj). # Check that assay is either RNA, SCT or integrated, not CSP

# Run Azimuth - run in RStudio
# - if a path to an Azimuth reference dir in not provided, use:
#    reference = "pbmcref"
# - to hard code the path to a local reference use:
#    reference = ref_dir_azimuth
results_azimuth  <- RunAzimuth(seurat.obj, reference = "pbmcref")

# Add annotations back into Seurat object
celltypes_azimuth <- results_azimuth@meta.data
annot_cols        <- names(celltypes_azimuth)
annot_cols        <- annot_cols[annot_cols %in% c("predicted.celltype.l1","predicted.celltype.l2","predicted.celltype.l3")]
celltypes_azimuth <- celltypes_azimuth[,annot_cols]

seurat.obj@meta.data <- cbind(seurat.obj@meta.data,celltypes_azimuth)


# 3. UMAP celltype annotations ----

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


# SAVE ----
saveRDS(seurat.obj,file = fh_processed_seurat_obj)
