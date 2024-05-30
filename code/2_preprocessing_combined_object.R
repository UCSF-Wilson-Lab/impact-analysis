#!/usr/bin/env Rscript

# START ----
# Workflow for processing combined analysis object
#  - Prior to running this code set working directory to 'impact-analysis/code'

library(Seurat)
library(rjson)
library(tidyverse)
library(stringr)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(SingleR)
library(kableExtra)
library(devtools)
library(parallel)
library(data.table)
library(harmony)
library(DoubletFinder)
library(impactSingleCellToolkit)

# INPUT and OUTPUT Directories
param_file_fh = "../input/input_one_patient_analysis.json"
params        = fromJSON(file = param_file_fh)

# INPUT
metadata_fh   = params$metadata_file
metadata      = read.csv(metadata_fh,stringsAsFactors = F,check.names = F)
THREADS       = params$threads

# OUTPUT Results Directories
plot_dir     = params$plot_directory
objects_dir  = params$objects_directory
if(!file.exists(plot_dir)){dir.create(plot_dir,recursive = TRUE)}
if(!file.exists(objects_dir)){dir.create(objects_dir,recursive = TRUE)}

# Input and Output RData Objects
fh_raw_seurat_obj       <- file.path(objects_dir,"seurat_obj.combined.RData")
fh_processed_seurat_obj <- file.path(objects_dir,"seurat_obj.processed.rds")

# Increase size of ENV
options(future.globals.maxSize= 891289600)
set.seed(112358)



# 1. Load and start normalization/clustering ----

# Load combined GEX/CSP object
load(fh_raw_seurat_obj)

# Omit Doublets 
cells_to_keep  <- names(seurat.obj$doublet_finder[seurat.obj$doublet_finder %in% "Singlet"])
seurat.obj     <- seurat.obj[,colnames(seurat.obj) %in% cells_to_keep]

# Re-cluster
seurat.obj   <- scaleAndClusterSeuratObject(seurat.obj,dims = 1:30,npca = 10,tsne = T)
pc_elbowplot <- plotOptimalPCsforSeuratObject(seurat.obj)



# SAVE ----
saveRDS(seurat.obj,file = fh_processed_seurat_obj)
