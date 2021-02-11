library(cluster)
library(parallel)
library(parallelDist)


calculate_silhouette <- function(dist.m, cell_identity, summary_stat = median, nThreads = 1){
    cell_cluster_int = as.numeric(cell_identity)
    b = closest_dist_other_cells = unlist(mclapply(1:nrow(dist.m), 
                                          function(i) {
                                              this_cell = which(cell_identity == cell_identity[i])
                                              stopifnot(length(this_cell) == 2)
                                              dist.m[i, which.min(dist.m[i, -this_cell])]
                                          }
                                          , mc.cores = nThreads))
    a = dist_self = get_self_reflection_neighbor_rank(dist.m, cell_identity)
    
    #res = silhouette(cell_cluster_int, dist.m)
    sidx = (b - a) / pmax(a, b)
    return(summary_stat(sidx))
}

get_self_reflection_neighbor_rank <- function(rank_matrix, cell_identity){
    self_reflection_neighbor_rank = unlist(sapply(1:nrow(rank_matrix), 
        function(i) {
            cellids = rank_matrix[i, which(cell_identity == cell_identity[i])]; 
            ranks = cellids[cellids > 0]
            stopifnot(length(ranks) == 1)
            ranks[1]
        }

        ))
    stopifnot(length(self_reflection_neighbor_rank) == length(cell_identity))
    return(self_reflection_neighbor_rank)
}

calculate_self_dist_rank <- function(dist.rank.m, cell_identity, summary_stat = median){
    ks = get_self_reflection_neighbor_rank(dist.rank.m, cell_identity)
    k_frac = ks / length(ks)
    return(summary_stat(ks))
}

calculate_anchoring_dist <- function(dist.m, cell_identity, summary_stat = median){
    anchoring_dist = get_self_reflection_neighbor_rank(dist.m, cell_identity)
    cells = colnames(dist.m)
    modality = sapply(strsplit(cells, split = ".", fixed =TRUE), function(s) s[1])
    stopifnot(modality %in% c("atac", "rna"))
    # Normalize by the max distance in other modality
    dist.other_modality.m = dist.m
    for (mod in unique(modality)){
        dist.other_modality.m[modality == mod, modality == mod] = 0
    }
    norm_factor = apply(dist.other_modality.m, 1, max)
    anchoring_dist = anchoring_dist / norm_factor
    return(summary_stat(anchoring_dist))
}

calculate_ARI <- function(embedding, k = 50){
    ks = kmeans(embedding[,1:nDims], centers = k, nstart = 1)
   in_same_celltype = sapply(levels(embedding$cell.ft), function(lvl){
                          cell_types = ks$cluster[which(embedding$cell.ft == lvl)]
                          if( length(cell_types) != 2) {print(lvl); stop()}
                          return( as.numeric(cell_types[1] == cell_types[2]) )})
#    res = mclust::adjustedRandIndex(, embedding$cell.ft)
   print(str(in_same_celltype))
    return(mean(in_same_celltype))
}

calculate_graph_connectivity <- function(dist.rank.m, cell_identity, k = 50){
    knn_graph = apply(dist.rank.m, 1, function(x) ifelse(x <= k & x > 1, 1, 0))
    ks = get_self_reflection_neighbor_rank(dist.rank.m, cell_identity)
    return(mean(ks <= k))
}

get_summarisable_metrics <- function(dist.m, dist.rank.m, embedding, cell_identity, summary_stat){
    summary_sil = calculate_silhouette(dist.m, embedding$cell.ft, summary_stat = summary_stat)
    summary_dist_rank = calculate_self_dist_rank(dist.rank.m, embedding$cell.ft, summary_stat = summary_stat)
    summary_anchor_dist = calculate_anchoring_dist(dist.m, embedding$cell.ft, summary_stat = summary_stat)
    res.df = data.frame(Silhouette = summary_sil, 
                        Anchoring_dist_rank = summary_dist_rank, 
                        Anchoring_dist = summary_anchor_dist, 
                        summary_stat = as.character(substitute(summary_stat)))
    message("Done calculating summarisable metrics")
    print(res.df)
    return(res.df)  
}

get_value_metrics <- function(dist.m, dist.rank.m, embedding, cell_identity, summary_stat){
    metrics = get_summarisable_metrics(dist.m, dist.rank.m, embedding, cell_identity, summary_stat)
    metrics$Graph_connectivity = calculate_graph_connectivity(dist.rank.m, embedding$cell.ft, k = 50)
    metrics$ARI = calculate_ARI(embedding)
    metrics$summary_method = as.character(substitute(summary_stat))
    return(metrics)
} 

write_value_metrics <- function(dist.m, dist.rank.m, embedding, cell_identity){
    medians = get_value_metrics(dist.m, dist.rank.m, embedding, cell_identity, median)
    means = get_value_metrics(dist.m, dist.rank.m, embedding, cell_identity, mean)
    write.table(rbind(medians, means), file = paste0(output_path,"/summary_metrics.txt"), 
                quote = F, row.names = T, col.names= T, sep = "\t")
    message("wrote summary metrics")
}

write_cell_metrics <- function(dist.m, dist.rank.m, embedding, cell_identity){
    res.df = get_summarisable_metrics(dist.m, dist.rank.m, embedding, cell_identity, identity)
    write.table(res.df, file = paste0(output_path,"/per_cell_metrics.txt"), 
                quote = F, row.names = T, col.names = T, sep = "\t")
    message("wrote per cell metrics")
}

