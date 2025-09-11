#!/usr/bin/env Rscript

# START ----
# Code to plot hamming distance disctribution across immcantation results
#  - Prior to running this code set working directory to 'impact-analysis/code/0_vdj_alignments'

# Immcantation
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(shazam))
suppressPackageStartupMessages(library(tigger))
suppressPackageStartupMessages(library(dowser))
suppressPackageStartupMessages(library(airr))

suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(glmGamPoi))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dittoSeq))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(impactSingleCellToolkit))

# Increase size of ENV
options(future.globals.maxSize= 40000*1024^2)
set.seed(12345)

# 1. Load inputs

# INPUT

# Immcantation Results directory
imm_results_bcr_dir <- "../../../Results/results_immcantation_full_cohort/BCR_CSFPB"
imm_results_tcr_dir <- "../../../Results/results_immcantation_full_cohort/TCR_CSFPB"

# OUTPUT Results Directories
table_dir  = "../../../Results/results_immcantation_full_cohort"

# 2. Load and format TCR/BCR repertoire results
# * Clonotype definitions for BCR are based on the Heavy chain only
# * Plots CDR3 distance distribution among BCRs and come up with a new threshold
# * Clonotype definitions for TCR include both Alpha and Beta chains

### a. Load Immcantation Results

# File paths
bcr_heavy_fh <- file.path(imm_results_bcr_dir,"BCR_CSFPB_heavy_germ-pass.tsv")
bcr_light_fh <- file.path(imm_results_bcr_dir,"BCR_CSFPB_light_germ-pass.tsv")
tcr_beta_fh  <- file.path(imm_results_tcr_dir,"TCR_CSFPB_heavy_germ-pass.tsv")
tcr_alpha_fh <- file.path(imm_results_tcr_dir,"TCR_CSFPB_light_germ-pass.tsv")

# TCR
db.tra   <- read_rearrangement(tcr_alpha_fh) %>% as.data.frame()
db.trb   <- read_rearrangement(tcr_beta_fh) %>% as.data.frame()
# BCR
db.heavy <- read_rearrangement(bcr_heavy_fh) %>% as.data.frame()
db.light <- read_rearrangement(bcr_light_fh) %>% as.data.frame()


### b. filter contigs and merge all chains
# * All contigs are in frame
bcr_df <- mergeVDJresults(db.heavy,db.light,umi.thresh = 2,assay = "bcr")
tcr_df <- mergeVDJresults(db.trb,db.tra,umi.thresh = 2, assay = "tcr")

# format VDJ unique cell IDs to match scRNA-Seq (function only for this workshop)
formatUniqueIDcolum <- function(uni_id_column,sep = ":") {
  samples <- tstrsplit(uni_id_column, sep)[[1]]
  cells   <- tstrsplit(uni_id_column, sep)[[2]]
  fmt_ids <- paste(cells,samples,sep = "-")
  
  return(fmt_ids)
}

bcr_df$unique_cell_id <- formatUniqueIDcolum(bcr_df$unique_cell_id)
tcr_df$unique_cell_id <- formatUniqueIDcolum(tcr_df$unique_cell_id)


# Write TSVs
write.table(bcr_df,file = file.path(imm_results_bcr_dir,"IMPACT_ALL_BCR_immcantation_results.tsv"),row.names = F,sep = "\t",quote = F)
write.table(tcr_df,file = file.path(imm_results_tcr_dir,"IMPACT_ALL_TCR_immcantation_results.tsv"),row.names = F,sep = "\t",quote = F)

