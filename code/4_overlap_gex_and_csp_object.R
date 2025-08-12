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

# INPUT and OUTPUT Directories
# 4_overlap_gex_and_csp_object.R  [Param File]

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

# OUTPUT Results Directories
plot_dir     = params$plot_directory
objects_dir  = params$objects_directory
if(!file.exists(plot_dir)){dir.create(plot_dir,recursive = TRUE)}
if(!file.exists(objects_dir)){dir.create(objects_dir,recursive = TRUE)}

# Input Object
fh_processed_seurat_obj <- file.path(objects_dir,"seurat_obj.processed.rds")

# RData Objects to save
fh_overlapped_seurat_obj <- file.path(objects_dir,"seurat_obj_csp_overlap.processed.rds")

# Increase size of ENV
options(future.globals.maxSize= 40000*1024^2)


# 1. Separate samples by batch/assay ----

# create batch:[gex/csf/multi] column
metadata$run_group <- paste(metadata$run,metadata$type,sep = ":")
metadata_gex       <- metadata[metadata$type %in% c("counts","counts_gex"),]
metadata_csp       <- metadata[metadata$type %in% c("counts","counts_csp"),]

run_metadata_csp_list <- split(metadata_csp,metadata_csp$run_group)


# 2. Create Seurat Objects for each batch ----

# makeRunInputMtx
# - output.type = c("gex","csp")
makeRunInputMtx <- function(runID ,run_metadata_list,THREADS,
                            output.type = "gex",min.genes.gex=700,min.genes.csp=9) {
  run          <- tstrsplit(runID,":")[[1]]
  type         <- tstrsplit(runID,":")[[2]]
  multi.status <- FALSE
  if(type == "counts"){
    multi.status <- TRUE
  }
  
  run_metadata <- run_metadata_list[[runID]]
  run_metadata <- run_metadata[run_metadata$type %in% type,]
  
  # Get run directory and sample vector
  sample_dir_vec <- run_metadata$results_directory_path
  dataset_loc    <- tstrsplit(sample_dir_vec,"multi_counts/")[[1]] %>% unique()
  dataset_loc    <- paste(dataset_loc,"multi_counts/",sep = "")
  
  # Separate samples based on whether they are mult, gex only or csp only
  samples.vec   <- run_metadata$sample
  
  if(THREADS > length(samples.vec)){
    THREADS <- length(samples.vec)
  }
  min_genes <- 0
  if(output.type == "gex"){min_genes <- min.genes.gex}
  if(output.type == "csp"){min_genes <- min.genes.csp}
  gex.matrix <- generateCombinedMatrix(dataset_loc, samples.vec,THREADS = THREADS,multi.results = multi.status,
                                       assay = output.type,min.genes.per.cell = min_genes,max.genes.per.cell = NULL)
  
  return(gex.matrix)
}

csp_run_list <- names(run_metadata_csp_list) %>% as.list()

# Initialize variables
dataset_loc        <- ""
samples.vec        <- c()
multi.results      <- NULL
assay              <- NULL
min.genes.per.cell <- NULL
max.genes.per.cell <- NULL

# CSP input matrix
csp_mtx_list <- lapply(csp_run_list, makeRunInputMtx,
                       run_metadata_list=run_metadata_csp_list,
                       THREADS=THREADS,output.type = "csp",min.genes.csp=1)
merged_csp_mtx <- do.call("cbind",csp_mtx_list)


# 3. Load Processed object ----
seurat.obj <- readRDS(fh_processed_seurat_obj)


# 4. Overlap objects ----
cells_csp <- colnames(merged_csp_mtx)
cells_gex <- colnames(seurat.obj)
overlap   <- cells_csp[cells_csp %in% cells_gex]

# Only keep overlapping cells
seurat.obj <- seurat.obj[,colnames(seurat.obj) %in% overlap]
counts_csp <- merged_csp_mtx[,colnames(merged_csp_mtx) %in% overlap]

seurat.obj[["CSP"]] <- CreateAssayObject(counts = counts_csp)


# SAVE ----
saveRDS(seurat.obj,file = fh_overlapped_seurat_obj)

