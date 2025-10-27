#!/usr/bin/env Rscript

# START ----
# Workflow for processing combined analysis object
#  - Prior to running this code set working directory to 'impact-analysis/code/batched_celltype_annot'

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
# 3_split_combined_object_for_celltype_annot.R  [Param File]

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("ERROR: At least one argument must be supplied (JSON parameter file).json", call.=FALSE)
} 
param_file_fh = args[1]

#param_file_fh = "../../input/input_all_csf_analysis.json"
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

# 2. Create column to batch the object (3 object subsets) ----

# figure out the grouping based on run size
run_summary <- table(seurat.obj$batch) %>% as.data.frame() %>% setNames(c("run","cell_count"))
run_summary <- run_summary[order(run_summary$cell_count,decreasing = T),]
run_summary$object_subset        <- "object5"
run_summary[1,"object_subset"]   <- "object1" # largest batch
run_summary[2,"object_subset"]   <- "object2"
run_summary[3:4,"object_subset"] <- "object3"
run_summary[5:7,"object_subset"] <- "object4"

runs_subset1 <- run_summary[run_summary$object_subset %in% "object1","run"] %>% as.character()
runs_subset2 <- run_summary[run_summary$object_subset %in% "object2","run"] %>% as.character()
runs_subset3 <- run_summary[run_summary$object_subset %in% "object3","run"] %>% as.character()
runs_subset4 <- run_summary[run_summary$object_subset %in% "object4","run"] %>% as.character()
runs_subset5 <- run_summary[run_summary$object_subset %in% "object5","run"] %>% as.character()

# Create column in the big seurat object
seurat.obj$object_subset <- seurat.obj$batch
seurat.obj$object_subset[seurat.obj$object_subset %in% runs_subset1] <- "object1"
seurat.obj$object_subset[seurat.obj$object_subset %in% runs_subset2] <- "object2"
seurat.obj$object_subset[seurat.obj$object_subset %in% runs_subset3] <- "object3"
seurat.obj$object_subset[seurat.obj$object_subset %in% runs_subset4] <- "object4"
seurat.obj$object_subset[seurat.obj$object_subset %in% runs_subset5] <- "object5"


# 2. Split and save objects ----
seurat.obj.list <- SplitObject(seurat.obj, split.by="object_subset")

for (i in 1:length(seurat.obj.list)) {
  file_name <- paste("TEMP_object",i,".rds",sep = "")
  fh_subset_processed_seurat_obj <- file.path(objects_dir,file_name)
  
  ### SAVE RDS
  saveRDS(seurat.obj.list[[i]],file = fh_subset_processed_seurat_obj)
}
