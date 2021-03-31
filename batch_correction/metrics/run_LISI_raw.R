library(ggplot2)
library(lisi)
library(ggpubr)

# Script for visualization 
rm(list=ls())
args = commandArgs(trailingOnly = TRUE)

source('metrics/lisi_utils.R')

method_use = args[1]
pca_file = args[2]
output_dir = args[3]
nPCs = as.integer(args[4])

#this_dir <- '/home/hoa/hoatran/demo_normalization/manuscript_results/dataset2_cellatlas/'
#setwd(this_dir)
eval_metric <- 'LISI' 
#dataset_use <- 
plx = 40

# Get output of LISI
#data_dir = paste0('/acrc/jinmiao/CJM_lab/hoatran/demo_normalization/xiaomeng/generate_PCA_tSNE_UMAP/',dataset_use,'/')

#print(fn_ls)
#dir.create(paste0(this_dir, eval_metric), showWarnings = F)


run_LISI_final('', pca_file, output_dir, eval_metric, method_use, plx, nPCs = nPCs)

#####################################
# iLISI batch
# Combine all results together
# 13 methods (except BBKNN) + raw data
#####################################
# View(head(resnet_df))
# # Important: With resnet, need to check the cell name 
# resnet_df$cell <- gsub('-[0-9]$','',resnet_df$cell)
# rownames(resnet_df) <- resnet_df$cell


#meta_ls <- get_celltype_common(pca_file)
#length(meta_ls$cells_common)

#get_cells_integration_iLISI_v2(dataset_use, meta_ls, this_dir, plx, eval_metric)



# #######################################
######### cLISI
# #######################################

#get_celltype_mixing_cLISI(dataset_use, this_dir, plx, eval_metric)

############################
## Summary output
############################


#median_LISI(output_dir, methods_use = method_use, plx = plx, eval_metric = eval_metric, nPCs = nPCs)
summary_LISI(output_dir, methods_use = method_use, plx = plx, eval_metric = eval_metric, nPCs = nPCs)

# fn <- 'dataset2_raw_pca.csv'
# meta_ls <- get_celltype_common(data_dir, fn)
# length(meta_ls$ct_common)


# summary_LISI <- function(meta_ls, this_dir, plottitle='LISI - dataset', plx=40, eval_metric='LISI/',ht=400, wd=400){
#   iLISI_df <- read.csv(paste0(this_dir, eval_metric,"result/",plx,"/iLISI_summary.csv"), head=T, check.names = F)
#   colns <- c('methods_use', 'iLISI_median', 'iLISI_median_norm')
#   
#   median_iLISI <- normalize_values(iLISI_df, colns, min_val=min(iLISI_df), max_val=max(iLISI_df))
#   mini <- min(median_iLISI$iLISI_median_norm)
#   maxi <- max(median_iLISI$iLISI_median_norm)
#   median_iLISI$iLISI_median_norm2 <- (median_iLISI$iLISI_median_norm - mini) / (maxi - mini)
#   
#   
#   cLISI_df <- read.csv(paste0(this_dir,eval_metric,"result/",plx,"/cLISI_summary.csv"), head=T, check.names = F)
#   colns <- c('methods_use', 'cLISI_median', 'cLISI_median_norm')
#   median_cLISI <- normalize_values(cLISI_df, colns, min_val=min(cLISI_df), max_val=max(cLISI_df))
#   median_cLISI$cLISI_median_norm_sub <- 1 - median_cLISI$cLISI_median_norm
#   minc <- min(median_cLISI$cLISI_median_norm_sub)
#   maxc <- max(median_cLISI$cLISI_median_norm_sub)
#   median_cLISI$cLISI_median_norm_sub2 <- (median_cLISI$cLISI_median_norm_sub - minc) / (maxc - minc)
#   
#   
#   final_df = merge(median_iLISI, median_cLISI, by="methods_use")
#   final_df$sum_normXY <- final_df$iLISI_median_norm2 + final_df$cLISI_median_norm_sub2
#   final_df$fscore <- (2 * final_df$iLISI_median_norm2 * final_df$cLISI_median_norm_sub2)/
#                      (final_df$iLISI_median_norm2 + final_df$cLISI_median_norm_sub2)
#   
#   final_df <- final_df[order(final_df$fscore, decreasing = T),]
#   # final_df$cLISI_median_norm_sub = 1 - final_df$cLISI_median_norm
#   write.csv(final_df,paste0(this_dir, eval_metric, "result/", plx, '/', 'summary_median_', plx, '.csv'), 
#               quote=F, row.names=F)
#   
# 
#   # plot final LISI
#   plot_final_LISI(final_df, plottitle, this_dir, plx, ht, wd, eval_metric, 
#                   xstring = 'cLISI_median_norm_sub', ystring = 'iLISI_median_norm', plottype = 'methods_use')
#   
# }

