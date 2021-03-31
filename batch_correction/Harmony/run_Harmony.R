
# Author : Kok Siong Ang 
# Date : 03/09/2019
# Proj : Harmony pipeline

########################
#load packages

library(harmony)
library(cowplot)
library(Seurat)  # Seurat 2 version
library(magrittr)

rm(list=ls())

args = commandArgs(trailingOnly = TRUE)
if ( length(args) != 4) {
    stop("Usage: Rscript run_script.R expr_filename metadata_filename dataset_name save_dir")
}


########################
#settings

filter_genes = F
filter_cells = F
normData = T
Datascaling = T
regressUMI = F
min_cells = 10
min_genes = 300
norm_method = "LogNormalize"
scale_factor = 10000
b_x_low_cutoff = 0.0125
b_x_high_cutoff = 3
b_y_cutoff = 0.5
numVG = 300
npcs = 20
visualize = T
outfile_prefix = args[3]
save_obj = T

src_dir = "Harmony/"
working_dir = args[4]
dir.create(working_dir, recursive = TRUE)

expr_filename = args[1]
metadata_filename = args[2]

batch_label = "batchlb"
celltype_label = "CellType"

########################
# read data 

expr_mat <- read.table(file = gzfile(expr_filename),sep="\t",header=T,row.names=1,check.names = F)
metadata <- read.table(file = gzfile(metadata_filename),sep="\t",header=T,row.names=1,check.names = F)

colnames(metadata)[colnames(metadata) == 'ct' | colnames(metadata) == 'celltype'] <- 'CellType'

expr_mat <- expr_mat[, rownames(metadata)]

########################
# run pipeline

source(paste0(src_dir,'call_harmony.R'))
#setwd(working_dir)

b_seurat = harmony_preprocess(
                expr_mat, metadata, 
                filter_genes = filter_genes, filter_cells = filter_cells,
                normData = normData, Datascaling = Datascaling, regressUMI = regressUMI, 
                min_cells = min_cells, min_genes = min_genes, 
                norm_method = norm_method, scale_factor = scale_factor, 
                b_x_low_cutoff = b_x_low_cutoff, b_x_high_cutoff = b_x_high_cutoff, b_y_cutoff = b_y_cutoff, 
                numVG = numVG, npcs = npcs, 
                batch_label = batch_label, celltype_label = celltype_label)

b_seurat = call_harmony(b_seurat, batch_label, celltype_label, npcs, plotout_dir = working_dir, saveout_dir = working_dir, outfilename_prefix = outfile_prefix, visualize = visualize, save_obj = save_obj)




