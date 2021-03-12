
# Author : Kok Siong Ang 
# Date : 03/09/2019
# Proj : Run Seurat 3 pipeline

########################
#load packages

library(Seurat)  # Seurat 2 version
packageVersion('Seurat')
library(magrittr)
library(cowplot)

rm(list=ls())

args = commandArgs(trailingOnly = TRUE)
if ( length(args) != 4) {
    stop("Usage: Rscript run_script.R expr_filename metadata_filename dataset_name save_dir")
}

########################
#settings

normData = T
Datascaling = T
regressUMI = F
min_cells = 0 #10
min_genes = 0 #300
norm_method = "LogNormalize"
scale_factor = 10000
numVG = 300
nhvg = 2000
npcs = 20
visualize = T
outfile_prefix = args[3]
save_obj = T

src_dir = "Seurat3/"
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

source(paste0(src_dir,'call_seurat_3.R'))
#setwd(working_dir)

batch_list = seurat3_preprocess(
                expr_mat, metadata, 
                normData = normData, Datascaling = Datascaling, regressUMI = regressUMI, 
                min_cells = min_cells, min_genes = min_genes, 
                norm_method = norm_method, scale_factor = scale_factor, 
                numVG = numVG, nhvg = nhvg, 
                batch_label = batch_label, celltype_label = celltype_label)

dir.create("Raw")
saveRDS(batch_list, "Raw/batch_list.RDS")

batches = call_seurat3(batch_list, batch_label, celltype_label, npcs, plotout_dir = working_dir, saveout_dir = working_dir, outfilename_prefix = outfile_prefix, visualize = visualize, save_obj = save_obj)




