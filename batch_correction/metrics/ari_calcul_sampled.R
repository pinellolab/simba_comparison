# Author : Nicole Lee
# Date : 29/08/2019
# Purpose: Second function to be called in ARI pipeline
#          Performs subsampling, extracts common cells for batch score
#          calculation, calculates the ARI scores for batch and cell type
#          Returns individual text files for each batch-correction method

####################################
# ARI Adjusted Rand Index (sampled)
####################################
# Input: myData: 20 PCs and batch, celltype columns
# cpcs: vector containing column number of PCs in myData

cosine_dist = function(embedding){
    Matrix <- as.matrix(embedding)
    sim <- Matrix / sqrt(rowSums(Matrix * Matrix))
    sim <- sim %*% t(sim)
    D_sim <- as.dist(1 - sim)
    return(D_sim)
}

get_Louvain_clusters <- function(embedding, nClusters, max_steps = 100, dissim = "euclidean"){
    require(Seurat)
    if (dissim == "euclidean"){
        nn_graph = FindNeighbors(embedding)
    } else if (dissim == "cosine"){
        dist_mat = cosine_dist(embedding)
        nn_graph = FindNeighbors(object = dist_mat)
    } else{
        stop("Not a valid metric")
    }
    steps = 0
    current_min = 0
    current_max = 3
    while(steps < max_steps){
        current_resolution = current_min + ((current_max - current_min)/2)
        clusters = FindClusters(nn_graph$snn, resolution = current_resolution, verbose = FALSE)[,1]
        current_n = length(unique(clusters))
        #print(paste0("got ", current_n, " at resolution ", current_resolution))
        if (current_n > nClusters) current_max = current_resolution
        else if (current_n < nClusters) current_min = current_resolution
        else{
            print(paste0("Got ", current_n, " at resolution ", current_resolution))
            return(clusters)}
        steps = steps + 1
    }
    print("Cannot find the number of clusters")
    print(paste0("Clustering solution from last iteration is used:", current_n, " at resolution ", current_resolution))
    return(clusters)
}

get_kmeans_clusters <- function(embedding, nClusters, max_steps = 100, dissim = "euclidean"){
    if (! dissim %in% c("euclidean", "Euclidean")){
        print(paste0("For kmeans clustering method ", dissim, " distance not implemented"))
        stop()
    }
    clustering_result <- kmeans(x = embedding, centers=nClusters, iter.max = max_steps)
    return(clustering_result$cluster)
}

get_clusters <- function(embedding, nClusters, max_steps = 100, dissim = "euclidean", clustering_method = c("Louvain", "kmeans")){
    fn = switch(clustering_method, 
           "Louvain" = get_Louvain_clusters, 
           "kmeans" = get_kmeans_clusters)
    return(fn(embedding, nClusters, max_steps, dissim))
}


ari_calcul_sampled <- function(myData, cpcs, 
                               method_use='resnet',
                               base_name='', nbiters=30, 
					           celltypelb='celltype', batchlb='batch', 
                              emb_type = emb_type, dissim = c('euclidean', 'cosine'), clustering_method = c("kmeans", "Louvain"))
{
  library(NbClust)
  library(mclust)
  set.seed(0)
  
  # get number of unique cell types
  nbct <- length(unique(myData[,celltypelb]))
  
  # get vector of unique cell types
  ce_types<-unique(myData[,celltypelb])
  
  # run function 20 times, each time extract 80% of data
  nbiters <- nbiters
  percent_extract <- 0.8
  
  it <- c()
  total_ari_batch <- c()
  total_ari_celltype <- c()
  
  # start loop for 20 times
  for(i in 1:nbiters) {
    
	# select cells for the subsampled dataset
    selectedcells<-vector()
    for (g in 1:nbct){
      cellpool<-which(myData[,celltypelb]==ce_types[g])
      ori_nbcells<-length(cellpool)
      cells_extract<-sample(cellpool, size=round(ori_nbcells*percent_extract), replace = F)
      selectedcells<-c(selectedcells, cells_extract)
    }
   
    selectedcells<-sort(selectedcells)
    
	# create the subsampled dataset
    myPCAExt <- myData[selectedcells,]
    
    ###############################
    # Clustering
    ###############################
    
    clustering_result <- get_clusters(embedding = myPCAExt[,cpcs], nClusters = nbct, dissim = dissim, clustering_method = clustering_method)
    myPCAExt$clusterlb <- clustering_result
      
    
	# assign the current myPCAExt to a unique object so that it can be stored later 
    assign(paste0("myPCAExt",i), myPCAExt)
	
    # Following clustering, get list of common cell types
    mySample <- subset(myPCAExt,select=c(celltypelb, batchlb))
    batches <- unique(mySample[,batchlb])
    
    ctls <- list()
    count <- 0
    for (b in batches){
      count <- count + 1
      ct <- unique(mySample[which(mySample[,batchlb]==b), celltypelb])
      ctls[[count]] <- ct
    }
    
    for(t in rep(1:length(ctls))){
      if(t==1){
        ct_common <- intersect(ctls[[t]], ctls[[t+1]])    
      }
      if(t>2){   #more than 2 batches
        ct_common <- intersect(ct_common, ctls[[t]])    
      }
    }
	# ct_common: common cell types amongst all batches
    
    cells_common <- rownames(mySample)[which(mySample[,celltypelb] %in% ct_common)]
    
    # create dataset with only common cells
    smallData <- myPCAExt[cells_common,]
    #assign(paste0("smallData",i), smallData)
    
    ###############################
    # ARI
    ###############################
    
    # run ARI
    ari_batch <- mclust::adjustedRandIndex(smallData[,batchlb], smallData$clusterlb)
    ari_celltype<-mclust::adjustedRandIndex(myPCAExt[,celltypelb], myPCAExt$clusterlb)
    
    it <- c(it,i)
    total_ari_batch <- c(total_ari_batch, ari_batch)
    total_ari_celltype <- c(total_ari_celltype, ari_celltype)
  }   # End of loop 
  
  
  # once looped 20 times to produce a total of 40 scores, next step follows:
  it <- c(it,nbiters+1)
  total_ari_batch <- c(total_ari_batch, median(total_ari_batch))
  total_ari_celltype <- c(total_ari_celltype, median(total_ari_celltype))
  
  methods <- rep(method_use, nbiters)
  methods <- c(methods,paste0(method_use,'_median'))
  
  # create final dataframe containing raw and median ARI scores
  myARI <- data.frame("use_case"=methods, 
                      "iteration"=it,
                      "ari_batch"=total_ari_batch, 
                      "ari_celltype"=total_ari_celltype)
  
  # write final dataframe to a text file
  write.table(myARI, file = paste0(base_name,method_use,"_ARI_", emb_type, "_", dissim, ".txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
  
  print('Save output in folder')
  print(base_name)
  
  return(list(myARI, myPCAExt1, myPCAExt2, myPCAExt3, myPCAExt4, myPCAExt5, myPCAExt6,
              myPCAExt7, myPCAExt8, myPCAExt9, myPCAExt10, myPCAExt11, 
              myPCAExt12,  myPCAExt13, myPCAExt14, myPCAExt15, myPCAExt16,  
              myPCAExt17, myPCAExt18, myPCAExt19, myPCAExt20))
}
