# Author : Nicole Lee
# Date : 29/08/2019
# Purpose: Primary script to run ARI pipeline

# clear workspace
rm(list=ls())

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 4){
    stop("Usage: Rscript run_ARISampled.R software pca_filename out_dir nPCs")
}
# source relevant functions
#source("Batch-effect-removal-benchmarking/Script/evaluation/ARI/ARI_utils/run_ARISampled.R")
source("batch_correction/metrics/ari_calcul_sampled.R")
#source("batch_correction/Batch-effect-removal-benchmarking/Script/evaluation/ARI/ARI_utils/conclude_ARISampled.R")

##############################
############ Set arguments

eval_metric <- 'ARISampled'

# read in files from this_dir
pca_filename = args[2]
method_use = args[1]

# send output to out_dir
out_dir = args[3]
nPCs = as.integer(args[4])



# Author : Nicole Lee
# Date : 29/08/2019
# Purpose: First function to be called in ARI pipeline
#          Reads in relevant dataset and grabs essential columns
#          Calls the following function 'ari_calcul_sampled'
#'
#' @param fn '_pca.csv' file

run_ARISampled <- function(pca_file, out_dir, eval_metric, method_use){

  thisData <- read.table(pca_file, sep = "\t", head = T, row.names = 1)

  # Get relevant columns
  colPCA <- grep('([Pp][Cc]_?)|(V)|(harmony.?)|(W)|(D)',colnames(thisData))
  colPCA <- colPCA[1:nPCs]
  str(thisData)
  colnames(thisData)[grep('[cC]ell_?[tT]ype',colnames(thisData))] <- 'celltype'
  colnames(thisData)[grep('([bB]atch)|(BATCH)|(batchlb)',colnames(thisData))] <- 'batch'

  #setwd(paste0(out_dir, '/', eval_metric, "_OP"))
  #temp<-ari_calcul_sampled(myData=thisData, cpcs=colPCA, isOptimal=TRUE,
  #                         method_use = method_use,
  #                         base_name=paste0(dataset_no, eval_metric, '_OP_'))
  setwd(out_dir)
  print(colnames(thisData))
  temp<-ari_calcul_sampled(myData=thisData, cpcs=colPCA, isOptimal=FALSE,
                           method_use = method_use,
                           base_name='')
  return(temp)
}


Rseurat3<-run_ARISampled(pca_filename, out_dir, eval_metric, method_use)

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
