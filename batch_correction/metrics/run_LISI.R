library(ggplot2)
library(lisi)

# Script for visualization 
rm(list=ls())
args = commandArgs(trailingOnly = TRUE)

source('metrics/lisi_utils.R')

pca_file = args[1]
output_dir = args[2]
emb_type = args[3]
dataset_use = args[4]
method_use = c(args[5])

plx = as.numeric(args[6])



# Get output of LISI
message(paste0("Calculating LISI for ", paste0(method_use, collapse = ", ")))

#print(fn_ls)
#dir.create(paste0(this_dir, eval_metric), showWarnings = F)
eval_metric = "LISI"
lisi.df = run_LISI_final('', pca_file, output_dir, eval_metric, method_use, plx, emb_type = emb_type)

shared_celltype_cells <- get_celltype_common(pca_file)$cells_common
cLISI_common = lisi.df$cell_type[which(lisi.df$cell %in% shared_celltype_cells)]

# Get summary metrics (median)
outfile = paste0(output_dir, "/", method_use, "_LISI_", plx, "_", emb_type, ".txt")
res.df = data.frame(methods_use = method_use, 
                    iLISI_median = median(lisi.df$batch),
                    cLISI_median = median(cLISI_common), 
                    plx = plx)
write.table(res.df, file = outfile, quote = F, col.names = T, row.names = F, sep = " ")
