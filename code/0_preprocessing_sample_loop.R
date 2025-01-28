#!/usr/bin/env Rscript

# START ----
# Workflow for creating analysis objects
#  - Prior to running this code set working directory to 'impact-analysis/code'

library(Seurat)
library(rjson)
library(tidyverse)
library(stringr)
library(knitr)
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

library(SoupX)
library(dsb)
library(DropletUtils)

# INPUT and OUTPUT Directories
# 0_preprocessing_sample_loop.R  [Param File]

#args = commandArgs(trailingOnly=TRUE)
#if (length(args)!=1) {
#  stop("ERROR: At least one argument must be supplied (JSON parameter file).json", call.=FALSE)
#} 
#param_file_fh = args[1]

param_file_fh = "../input/input_small_patient_subset_analysis.json"
params        = fromJSON(file = param_file_fh)

# INPUT
metadata_fh   = params$metadata_file
metadata      = read.csv(metadata_fh,stringsAsFactors = F,check.names = F)
THREADS       = params$threads

# Functions ----
# These functions could be incorporated into the R package

### runAmbientRNAfilterAndOuputFiles
runAmbientRNAfilterAndOuputFiles <- function(raw_hd5_fh,filt_hd5_fh,output_dir) {
  # Load and only focus on Gene Expression
  raw_mtx  <- Read10X_h5(raw_hd5_fh,use.names = T)
  raw_mtx  <- raw_mtx$`Gene Expression`
  filt_mtx <- Read10X_h5(filt_hd5_fh,use.names = T)
  filt_mtx <- filt_mtx$`Gene Expression`
  
  # make object
  seurat_obj   <- CreateSeuratObject(counts = filt_mtx)
  soup_channel <- SoupChannel(raw_mtx, filt_mtx)
  
  # process object
  seurat_obj    <- SCTransform(seurat_obj, verbose = F)
  seurat_obj    <- RunPCA(seurat_obj, verbose = F)
  seurat_obj    <- RunUMAP(seurat_obj, dims = 1:30, verbose = F)
  seurat_obj    <- FindNeighbors(seurat_obj, dims = 1:30, verbose = F)
  seurat_obj    <- FindClusters(seurat_obj, verbose = T)
  
  # Add clusters to channel
  meta_obj    <- seurat_obj@meta.data
  umap_obj    <- seurat_obj@reductions$umap@cell.embeddings
  soup_channel  <- setClusters(soup_channel, setNames(meta_obj$seurat_clusters, rownames(meta_obj)))
  soup_channel  <- setDR(soup_channel, umap_obj)
  
  # Calculate ambient RNA
  soup_channel  <- autoEstCont(soup_channel)
  # Convert to integers
  adj_matrix  <- adjustCounts(soup_channel, roundToInt = T)
  
  # Output
  DropletUtils:::write10xCounts(output_dir, adj_matrix, version = "3")
}


### preprocessCountsUsingMetadata ----
preprocessCountsUsingMetadata <- function(sample,metadata_counts,ambient.rna.filter = TRUE) {
  metadata_sample <- metadata_counts[metadata_counts$sample %in% sample,]
  counts_dir      <- metadata_sample$results_directory_path
  hd5_dir         <- metadata_sample$path_hd5
  counts_filt_dir <- metadata_sample$path_counts_filt
  
  # Directory to read in and apply filters and doublet removal
  input_dir <- counts_dir
  if(ambient.rna.filter){input_dir <- counts_filt_dir}
  
  # If performing ambient RNA filter output results to filtered counts folder
  if(ambient.rna.filter){
    raw_mtx_fh  <- file.path(hd5_dir,"raw_feature_bc_matrix.h5")
    filt_mtx_fh <- file.path(hd5_dir,"sample_filtered_feature_bc_matrix.h5")
    
    if(! dir.exists(counts_filt_dir)){
      runAmbientRNAfilterAndOuputFiles(raw_hd5_fh = raw_mtx_fh,
                                       filt_hd5_fh = filt_mtx_fh,
                                       output_dir = counts_filt_dir)
    }
  }
  
  # Gene count filter
  
  # Doublet removal
  
}


# 1. list samples ----
metadata_counts <- metadata[metadata$type %in% "counts",]
metadata_counts$path_hd5 <- str_replace_all(metadata_counts$results_directory_path,"\\/multi_counts\\/","/multi_counts_hd5/")
metadata_counts$path_counts_filt <- str_replace_all(metadata_counts$results_directory_path,"\\/multi_counts\\/","/multi_counts_filt/")


samples <- metadata_counts$sample

# 2. loop through samples and pre-process ----
lapply(as.list(samples), preprocessCountsUsingMetadata,
       metadata_counts=metadata_counts,ambient.rna.filter=TRUE)
