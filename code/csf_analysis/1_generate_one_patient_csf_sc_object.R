#!/usr/bin/env Rscript

# START ----
# Workflow for creating analysis objects
#  - Prior to running this code set working directory to 'impact-analysis/code/csf_analysis'

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
patient_label = args[1]
param_file_fh = "../../input/input_all_csf_analysis.json"
params        = fromJSON(file = param_file_fh)

# INPUT
metadata_fh   = params$metadata_file
metadata      = read.csv(metadata_fh,stringsAsFactors = F,check.names = F)
THREADS       = params$threads

# OUTPUT Results Directories
plot_dir       = params$plot_directory
objects_dir    = params$objects_directory
pt_objects_dir = file.path(objects_dir,"patient_objects")
if(!file.exists(plot_dir)){dir.create(plot_dir,recursive = TRUE)}
if(!file.exists(objects_dir)){dir.create(objects_dir,recursive = TRUE)}
if(!file.exists(pt_objects_dir)){dir.create(pt_objects_dir,recursive = TRUE)}

source("functions_one_patient_sc_object.R",local = T)

# RData Objects to save
fh_raw_seurat_obj <- file.path(objects_dir,"seurat_obj_CSF.combined.RData")

# Increase size of ENV
options(future.globals.maxSize= 891289600)


# 1. Separate samples by batch/assay ----

# create batch:[gex/csf/multi] column
metadata$run_group <- paste(metadata$run,metadata$type,sep = ":")
uni_patient_ids    <- unique(metadata$StudyID)
# Get this working on a subset
uni_patient_ids    <- uni_patient_ids[1:2] %>% as.list()

### Loop through each patient and generate an object

# Initialize variables
dataset_loc        <- ""
samples.vec        <- c()
multi.results      <- NULL
assay              <- NULL
min.genes.per.cell <- NULL
max.genes.per.cell <- NULL

lapply(uni_patient_ids,createOnePatientObject,metadata=metadata,THREADS=THREADS,pt_objects_dir=pt_objects_dir)

