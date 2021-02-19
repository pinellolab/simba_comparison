library(cluster)
library(parallel)

calculate_pdist = function(X, Y){ 
    stopifnot(dim(X) == dim(Y))
    apply(X - Y, 1, function(v) sqrt(sum(v^2)) ) }

get_reflection_idx = function(cell_identity){
    myself_idx = sapply(1:length(cell_identity), function(i){
                        cid = cell_identity[i]; 
                        indexes = 1:length(cell_identity); 
                        myself = indexes[cell_identity == cid ]; 
                        myself[myself != i] })
    return(myself_idx)
}

calculate_silhouette <- function(embedding, cell_identity, summary_stat = median, nThreads = 1){
    require(FNN)
    # Find 2 closest neighbors
    res = get.knn(embedding, k = 2)$nn.index
    
    # get distance to closest neighbor
    closest_other_cells = res[,1]
    reflection_idx = get_reflection_idx(cell_identity)
    idx_myself_is_closest = which(1:length(res[,1]) == reflection_idx)
    closest_other_cells[idx_myself_is_closest] = res[,2]
    b = calculate_pdist(embedding, embedding[closest_other_cells,])

    # get distance to myself
    my_reflection = embedding[get_reflection_idx(cell_identity),]
    a = dist_self = calculate_pdist(embedding, my_reflection)
    
    #res = silhouette(cell_cluster_int, dist.m)
    sidx = (b - a) / pmax(a, b)
    return(summary_stat(sidx))
}


sampled_self_rank <- function(embedding, cell_identity, size = 2000, seed = 1){
    set.seed(seed)
    cell_identity = factor(cell_identity)
    sampled_cid = sample(levels(cell_identity), size = size)
    sampled_cells = cell_identity[cell_identity %in% sampled_cid]
    reflection_idx = get_reflection_idx(sampled_cells)
    emb = embedding[cell_identity %in% sampled_cid, ]
    dist_matrix = as.matrix(dist(emb))
    rank_matrix = apply(dist_matrix, 1, rank) - 1
    #self_rank = mapply(function(row, i) row[i], as.data.frame(rank_matrix), reflection_idx)
    self_rank = sapply(1:nrow(rank_matrix), function(i) rank_matrix[i, reflection_idx[i]])
    return(self_rank)
}

calculate_self_dist_rank <- function(embedding, cell_identity, threads = 20){
    require(foreach)
    R = 20
    cl <- parallel::makeForkCluster(threads)
    doParallel::registerDoParallel(cl)
    # Samples 5% of the cell identity to check the neighborhood ranks.
    sample_size = min(2000, length(cell_identity) %/% R)
    self_ranks.l = foreach(r=1:R, .combine = 'c') %dopar% {
        list(sampled_self_rank(embedding, cell_identity, sample_size, seed = r) / sample_size)
    }
    
    parallel::stopCluster(cl)
    return(sapply(self_ranks.l, median))
    #return(sapply(self_ranks.l, median))
}

calculate_self_dist_rank_all <- function(embedding, cell_identity){
    require(foreach)
        # Samples 5% of the cell identity to check the neighborhood ranks.
    sample_size = length(cell_identity) %/% 2
    self_ranks.l = sampled_self_rank(embedding, cell_identity, sample_size, seed = 1) / sample_size

    return(self_ranks.l)
                #return(sapply(self_ranks.l, median))
}


calculate_anchoring_dist <- function(embedding, cell_identity, k = 50, summary_stat = median){
    reflection_idx = get_reflection_idx(cell_identity)
    my_reflection = embedding[reflection_idx,]
    dist_self = calculate_pdist(embedding, my_reflection)

    # Get normalization factor: mean distance to kNN (k = 50)
    #res = get.knn(embedding, k)
    #norm_factor = apply(res$nn.dist, 1, mean)

    # Normalize by randomly sampled distance
    dists = sapply(1:5000, function(i) {pair = sample(cell_identity, 2); dist(embedding[pair[1],], embedding[pair[2],])}
    norm_factor = mean(dists)

    anchoring_dist = dist_self / norm_factor
    return(summary_stat(anchoring_dist))
}

calculate_ARI <- function(embedding, cell_identity){
    require(Seurat)
    nn_graph = FindNeighbors(embedding) # k = 20 by default, use all dimensions
    clusters = FindClusters(nn_graph$snn)[,1] #res = 0.8 by default
    in_same_celltype = sapply(levels(cell_identity), function(lvl){
                          cell_types = clusters[which(cell_identity == lvl)]
                          if( length(cell_types) != 2) {print(lvl); stop()}
                          return( as.numeric(cell_types[1] == cell_types[2]) )})
#    res = mclust::adjustedRandIndex(, embedding$cell.ft)
    print(str(in_same_celltype))
    return(mean(in_same_celltype))
}

calculate_graph_connectivity <- function(embedding, cell_identity, k = 50){
    require(FNN)
    message(paste0("Calculating graph connectivity for k = ", k))
    knn_graph = get.knn(embedding, k = k)$nn.index
    message("Done obtaining knn graph for graph connectivity")
    reflection_idx = get_reflection_idx(cell_identity)
     
    reflection_in_knn = sapply(1:length(reflection_idx), function(i) {return(reflection_idx[i] %in% knn_graph[i,])})
    #reflection_in_knn = mapply(function(row, i) i %in% row, as.data.frame(knn_graph), reflection_idx)
    return(mean(reflection_in_knn))
}

get_summarisable_metrics <- function(embedding, cell_identity, summary_stat){
    summary_sil = calculate_silhouette(embedding, embedding$cell.ft, summary_stat = summary_stat)
    summary_dist_rank = calculate_self_dist_rank(embedding, embedding$cell.ft, summary_stat = summary_stat)
    summary_anchor_dist = calculate_anchoring_dist(embedding, embedding$cell.ft, summary_stat = summary_stat)
    res.df = data.frame(Silhouette = summary_sil, 
                        Anchoring_dist_rank = summary_dist_rank, 
                        Anchoring_dist = summary_anchor_dist, 
                        summary_stat = as.character(substitute(summary_stat)))
    message("Done calculating summarisable metrics")
    print(res.df)
    return(res.df)  
}

get_value_metrics <- function(embedding, cell_identity, summary_stat){
    metrics = get_summarisable_metrics(embedding, cell_identity, summary_stat)
    metrics$Graph_connectivity = calculate_graph_connectivity(embedding, cell_identity, k = 50)
    metrics$ARI = calculate_ARI(embedding)
    metrics$summary_method = as.character(substitute(summary_stat))
    return(metrics)
} 

write_value_metrics <- function(embedding, cell_identity){
    medians = get_value_metrics(embedding, cell_identity, median)
    means = get_value_metrics(embedding, cell_identity, mean)
    write.table(rbind(medians, means), file = paste0(output_path,"/summary_metrics.txt"), 
                quote = F, row.names = T, col.names= T, sep = "\t")
    message("wrote summary metrics")
}

write_cell_metrics <- function(embedding, cell_identity){
    res.df = get_summarisable_metrics(embedding, cell_identity, identity)
    write.table(res.df, file = paste0(output_path,"/per_cell_metrics.txt"), 
                quote = F, row.names = T, col.names = T, sep = "\t")
    message("wrote per cell metrics")
}

