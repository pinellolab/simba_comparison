
# Author : Marion Chevrier
# Date : 30/09/2019
# Proj : Run LIGER pipeline

########################
#load packages

#devtools::install_github('MacoskoLab/liger')
library(scales)
library(liger)
library(Matrix)
library(Rtsne)
library(ggplot2)

rm(list=ls())

args = commandArgs(trailingOnly = TRUE)
if ( length(args) != 5) {
    stop("Usage: Rscript run_script.R expr_filename metadata_filename dataset_name save_dir")
}


########################
#settings

var.thresh = 0.1
k = as.numeric(args[5])
nrep = 3
visualize = T
outfile_prefix = args[3]
save_obj = T

src_dir = "LIGER/"
working_dir = args[4]

expr_filename = args[1]
metadata_filename = args[2]

batch_label = "batchlb"
celltype_label = "CellType"

########################
# read data 

expr_mat <- read.table(file = gzfile(expr_filename),sep="\t",header=T,row.names=1,check.names = F)
metadata <- read.table(file = gzfile(metadata_filename),sep="\t",header=T,row.names=1,check.names = F)

colnames(metadata)[colnames(metadata) == 'ct' | colnames(metadata) == "celltype"] <- 'CellType'

expr_mat <- expr_mat[, rownames(metadata)]

########################
# run pipeline

source(paste0(src_dir,'call_liger.R'))

liger_obj <- liger_preprocess(expr_mat, metadata, 
                              var.thresh=var.thresh,
                              batch_label = batch_label)

call_liger(liger_obj, metadata, batch_label, celltype_label, k = k, nrep = nrep, 
           plotout_dir = working_dir, saveout_dir = working_dir, outfilename_prefix = outfile_prefix, visualize = visualize, save_obj = save_obj)
