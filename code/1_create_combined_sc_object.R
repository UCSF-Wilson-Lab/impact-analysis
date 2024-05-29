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

# RData Objects to save
fh_raw_seurat_obj <- file.path(objects_dir,"seurat_obj.combined.RData")

# Increase size of ENV
options(future.globals.maxSize= 891289600)


# 1. Separate samples by batch/assay ----

# create batch:[gex/csf/multi] column
metadata$run_group <- paste(metadata$run,metadata$type,sep = ":")
metadata_gex       <- metadata[metadata$type %in% c("counts","counts_gex"),]
metadata_csp       <- metadata[metadata$type %in% c("counts","counts_csp"),]


run_metadata_gex_list <- split(metadata_gex,metadata_gex$run_group)
run_metadata_csp_list <- split(metadata_csp,metadata_csp$run_group)


# 2. Create Seurat Objects for each batch

# makeRunInputMtx
# - output.type = c("gex","csp")
makeRunInputMtx <- function(runID ,run_metadata_list,THREADS,output.type = "gex") {
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
  gex.matrix <- generateCombinedMatrix(dataset_loc, samples.vec,THREADS = THREADS,multi.results = multi.status,
                                       assay = output.type,min.genes.per.cell = 700,max.genes.per.cell = NULL)
  
  return(gex.matrix)
}

gex_run_list <- names(run_metadata_gex_list) %>% as.list()
csp_run_list <- names(run_metadata_csp_list) %>% as.list()

csp_mtx_list <- lapply(csp_run_list, makeRunInputMtx,
                       run_metadata_list=run_metadata_csp_list,
                       THREADS=THREADS,output.type = "csp")



