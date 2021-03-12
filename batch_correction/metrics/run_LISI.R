library(ggplot2)
library(lisi)
library(ggpubr)

# Script for visualization 
rm(list=ls())
args = commandArgs(trailingOnly = TRUE)

source('metrics/lisi_utils.R')

pca_files = args[1:(length(args)-2)]
output_dir = args[length(args)-1]
emb_type = args[length(args)]

#this_dir <- '/home/hoa/hoatran/demo_normalization/manuscript_results/dataset2_cellatlas/'
#setwd(this_dir)
eval_metric <- 'LISI' 
tmp = strsplit(output_dir, split= "/")[[1]]
dataset_use <- tmp[length(tmp)]
plx = 40

# Get output of LISI
#data_dir = paste0('/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/xiaomeng/generate_PCA_tSNE_UMAP/',dataset_use,'/')

#methods_use = c("Raw_nopp", "Raw_PCA", "Seurat3", "Harmony", "LIGER", "SIMBA")
methods_use = sapply(strsplit(pca_files, split="/"), function(vec) vec[1])
message(paste0("Calculating LISI for ", paste0(methods_use, collapse = ", ")))

fn_ls <- pca_files

stopifnot(length(methods_use) == length(fn_ls))

#print(fn_ls)
#dir.create(paste0(this_dir, eval_metric), showWarnings = F)

for (i in rep(1:length(methods_use), 1)){
    print(methods_use[i])
    run_LISI_final('', pca_files[i], output_dir, eval_metric, methods_use[i], plx, emb_type = emb_type)
}

#####################################
# iLISI batch
# Combine all results together
# 13 methods (except BBKNN) + raw data
#####################################
# View(head(resnet_df))
# # Important: With resnet, need to check the cell name 
# resnet_df$cell <- gsub('-[0-9]$','',resnet_df$cell)
# rownames(resnet_df) <- resnet_df$cell

pca_file = paste0("Raw_PCA/output/", dataset_use,'_Raw_PCA_pca.txt')
meta_ls <- get_celltype_common(pca_file)
length(meta_ls$cells_common)
get_cells_integration_iLISI_v2(dataset_use, methods_use, meta_ls, output_dir, plx, eval_metric)



# #######################################
######### cLISI
# #######################################

get_celltype_mixing_cLISI(dataset_use, methods_use, output_dir, plx, eval_metric)

############################
## Summary output
############################


summary_LISI(output_dir, methods_use = methods_use, plx = plx, eval_metric = eval_metric, emb_type = emb_type)

