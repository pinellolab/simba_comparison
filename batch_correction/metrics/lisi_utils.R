normalize_values <- function(myData, colns, min_val=1, max_val=2){
  med_val <- c()
  med_val_norm <- c()
  methods_use <- c()
  for(cn in colnames(myData)){
      methods_use <- c(methods_use, cn)
      mi <- median(myData[,cn])
      med_val <- c(med_val, mi)
      med_val_norm <- c(med_val_norm, (mi - min_val) / (max_val - min_val))
  }
  myDataNorm <- data.frame('X1'=methods_use, 'X2'=med_val, 'X3'=med_val_norm)
  colnames(myDataNorm) <- colns
  return(myDataNorm)
}


get_celltype_common <- function(data_path){
  myData <- read.table(data_path, head=T, row.names = 1, check.names = FALSE, sep = "\t")
  colnames(myData)[grep('[cC]ell_?[tT]ype',colnames(myData))] <- 'cell_type'
  colnames(myData)[grep('([bB]atch)|(BATCH)',colnames(myData))] <- 'batch'
  batches <- unique(myData$batch)
  celltypels <- unique(myData$cell_type)
  #print(celltypels)
  if(length(batches)>1){
    ctls <- list()
    count <- 0
    for (b in batches){
      count <- count + 1
      ct <- unique(myData[which(myData$batch==b),'cell_type'])
      ctls[[count]] <- ct
    }
    for(i in rep(1:length(ctls))){
      if(i==2){
        ct_common <- intersect(ctls[[i-1]], ctls[[i]])    
      }
      if(i>2){   #more than 2 batches
        ct_common <- intersect(ct_common, ctls[[i]])    
      }
    }
    ct_common <- unique(ct_common)
    print(paste0("Common cell types: ", paste(ct_common, collapse = " ")))
    cells_common <- rownames(myData)[which(myData$cell_type %in% ct_common) ]
    return(list('ct_common'=ct_common, 'cells_common'=cells_common, 'batches'=batches, 'celltypels'=celltypels))
  } else{
    return(NULL)
  } 
  
}


compute_lisi_annoy = function (X, meta_data, label_colnames, perplexity = 30, nn_eps = 0, metric = metric)
{
    N <- nrow(meta_data)
    metric = switch(metric, "cosine" = "cosine", "Euclidean" = "euclidean", "euclidean" = "euclidean")
    dknn <- Seurat:::AnnoyNN(data = as.matrix(X), metric = metric, k = perplexity * 3)
    lisi_df <- data.frame(matrix(NA, N, length(label_colnames)))
    lisi_df <- Reduce(cbind, lapply(label_colnames, function(label_colname) {
        labels <- data.frame(meta_data)[, label_colname, drop = TRUE]
        if (any(is.na(labels))) {
            message("Cannot compute LISI on missing values")
            return(rep(NA, N))
        }
        else {
            dknn$nn.idx <- dknn$nn.idx[, 2:ncol(dknn$nn.idx)]
            dknn$nn.dists <- dknn$nn.dists[, 2:ncol(dknn$nn.dists)]
            labels <- as.integer(factor(labels)) - 1
            n_batches <- length(unique(labels))
            simpson <- lisi:::compute_simpson_index(t(dknn$nn.dists),
                t(dknn$nn.idx) - 1, labels, n_batches, perplexity)
            return(1/simpson)
        }
    }))
    lisi_df <- as.data.frame(lisi_df)
    colnames(lisi_df) <- label_colnames
    row.names(lisi_df) <- row.names(meta_data)
    return(lisi_df)
}

run_LISI_final <- function(fn, data_path, save_dir, eval_metric = "LISI", methods_use, plx=40, emb_type = emb_type, dissim = dissim){
  
  # myPCA <- read.csv(paste0(data_path, fn), head=T, row.names = 1, check.names = FALSE)
  myPCA <- read.table(data_path, head = T, row.names = 1, check.names = FALSE, sep = "\t")
  cpcs <- grep('([Pp][Cc]_?)|(D)|(V)|(harmony.?)|(W)|(UMAP)',colnames(myPCA))
  lisi_embeddings <- myPCA[,cpcs]
  
  colnames(myPCA)[grep('[cC]ell_?[tT]ype',colnames(myPCA))] <- 'cell_type'
  colnames(myPCA)[grep('([bB]atch)|(BATCH)|(batchlb)',colnames(myPCA))] <- 'batch'
  
  lisi_meta_data <- subset(myPCA, select=c('batch','cell_type'))
  
  lisi_label = c('batch', 'cell_type')
  
  lisi_res <- compute_lisi_annoy(lisi_embeddings, lisi_meta_data, lisi_label, perplexity = plx, metric = dissim)
  lisi_res$cell <- rownames(lisi_embeddings)
    
  shared_celltype_cells <- get_celltype_common(pca_file)$cells_common
  print(paste0(length(shared_celltype_cells), " shared cells"))
  lisi_res = lisi_res[which(lisi_res$cell %in% shared_celltype_cells),]

  lisi_batch <- subset(lisi_res,select=c('batch','cell'))
  lisi_celltype <- subset(lisi_res,select=c('cell_type','cell'))
  
  dir.create( paste0(save_dir, '/lisi_tmpdir/') )
  write.table(lisi_batch, paste0(save_dir, '/lisi_tmpdir/', methods_use, '_', eval_metric, '_batch_',plx,'.txt'), quote=F, sep='\t', row.names=T, col.names=NA)
  write.table(lisi_celltype,paste0(save_dir, '/lisi_tmpdir/', methods_use, '_', eval_metric, '_celltype_',plx,'.txt'), quote=F, sep='\t', row.names=T, col.names=NA)
  
  require(dplyr) 
  res_df = inner_join(lisi_batch, lisi_celltype, by = "cell")
  return(res_df)
}

