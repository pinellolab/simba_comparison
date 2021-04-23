# Author : Nicole Lee
# Date : 29/08/2019
# Purpose: Primary script to run ARI pipeline

# clear workspace
rm(list=ls())

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 4){
    stop("Usage: Rscript run_ARISampled.R software pca_filename out_dir emb_type similarity")
}
# source relevant functions
#source("Batch-effect-removal-benchmarking/Script/evaluation/ARI/ARI_utils/run_ARISampled.R")
source("metrics/ari_calcul_sampled.R")
#source("Batch-effect-removal-benchmarking/Script/evaluation/ARI/ARI_utils/conclude_ARISampled.R")

##############################
############ Set arguments

eval_metric <- 'ARISampled'

# read in files from this_dir
pca_filename = args[2]
method_use = args[1]

# send output to out_dir
out_dir = args[3]
emb_type = args[4]
dissim = args[5]


# Author : Nicole Lee
# Date : 29/08/2019
# Purpose: First function to be called in ARI pipeline
#          Reads in relevant dataset and grabs essential columns
#          Calls the following function 'ari_calcul_sampled'
#'
#' @param fn '_pca.csv' file

run_ARISampled <- function(pca_file, out_dir, eval_metric, method_use, nbiters,  emb_type = emb_type, dissim = dissim, clustering_method){

  thisData <- read.table(pca_file, sep = "\t", head = T, row.names = 1)

  # Get relevant columns
  colPCA <- grep('([Pp][Cc]_?)|(V)|(harmony.?)|(W)|(D)|(UMAP)',colnames(thisData))
  
  colnames(thisData)[grep('[cC]ell_?[tT]ype',colnames(thisData))] <- 'celltype'
  colnames(thisData)[grep('([bB]atch)|(BATCH)|(batchlb)',colnames(thisData))] <- 'batch'

  setwd(out_dir)
  temp<-ari_calcul_sampled(myData=thisData, cpcs=colPCA,
                           method_use = method_use,
                           nbiters = nbiters,
                           base_name='', emb_type = emb_type, dissim = dissim, 
                           clustering_method = clustering_method)
  return(temp)
}


res <-run_ARISampled(pca_filename, out_dir, eval_metric, method_use, nbiters = 20, emb_type = emb_type, dissim = dissim, clustering_method = "kmeans")

##############################
############ Extracting all data from all methods in dataset 

# Reads files from dir_this 
#dir_this<-paste0(out_dir, eval_metric, "_CT")
#dir_this<-paste0(out_dir, eval_metric, "_OP")

# Perform the following function to produce final CSV file
#wholedf<-conclude_ARISampled(dir_this, "Dataset_2") 

#save.image(paste0("Dataset2", "_complete.RData"))


##############################
############ Consolidate all raw data into one file

#setwd(dir_this)
#rm(list=ls())

#source("./ARI_utils/ARI_files_consolidate.R")
 
#ari_consolidate()

### END
