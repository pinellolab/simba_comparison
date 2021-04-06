library(cluster)
library(parallel)

calculate_pdist = function(X, Y){ 
    stopifnot(dim(X) == dim(Y))
    res = apply(X - Y, 1, function(v) sqrt(sum(v^2)) )
    return(res)
 }

get_reflection_idx = function(cell_identity){
    myself_idx = sapply(1:length(cell_identity), function(i){
                        cid = cell_identity[i]; 
                        indexes = 1:length(cell_identity); 
                        myself = indexes[cell_identity == cid ]; 
                        myself[myself != i] })
    return(myself_idx)
}


get_dotsim_dist_matrix <- function(embedding){
    embedding.m = as.matrix(embedding)
    sim_matrix = embedding.m %*% t(embedding.m)
    #dist_matrix = do.call(rbind, lapply(1:nrow(sim_matrix), function(i) {sim_matrix[i,i] - sim_matrix[i,]}))
    dist_matrix = do.call(rbind, lapply(1:nrow(sim_matrix), function(i) {max(sim_matrix[i,]) - sim_matrix[i,]}))
    dist_matrix[dist_matrix < 0] = 0
    return(dist_matrix)
}


calculate_silhouette <- function(embedding, cell_identity, summary_stat = median, nThreads = 1, 
                                 dissim = c("Euclidean", "dot product")){
    require(FNN)
    # Find 2 closest neighbors
    if (dissim == "Euclidean") {
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

    } else if (dissim == "dot product") {
        embedding.m = as.matrix(embedding)
        dots = embedding.m %*% t(embedding.m)
        self_sims = diag(dots)
        diag(dots) <- NA
        max_sims = apply(dots, 1, max, na.rm = TRUE)
        sidx = (self_sims - max_sims) / pmax(self_sims, max_sims)
    } else {
        stop("Invalid dissim in Silhouette")
    }
    
    return(summary_stat(sidx))
}


sampled_self_rank <- function(embedding, cell_identity, modality, size = 2000, seed = 1, threads = 20, dissim = c("Euclidean", "dot product")){
    set.seed(seed)
    require(parallel)
    cell_identity = factor(cell_identity)
    
    sampled_cid = sample(levels(cell_identity), size = size)
    sampled_cells = cell_identity[cell_identity %in% sampled_cid]
    
    reflection_idx = get_reflection_idx(sampled_cells)
    emb = embedding[cell_identity %in% sampled_cid, ]
    
    if (dissim == "Euclidean"){
    dist_matrix = as.matrix(dist(emb))
    rank_matrix = do.call(rbind, mclapply(1:nrow(dist_matrix), function(i) {
                                     this_mod = modality[i]
                                     dist_row = dist_matrix[i,]
                                     dist_row[which(modality == this_mod)] <- NA
                                     return(rank(dist_row) - 1)}, 
                                mc.cores = threads))
    } else{
        embedding.m = as.matrix(embedding)
        sim_matrix = embedding.m %*% t(embedding.m)
        rank_matrix = do.call(rbind, mclapply(1:nrow(sim_matrix), function(i) {
                                     this_mod = modality[i]
                                     sim_row = sim_matrix[i,]
                                     sim_row[which(modality == this_mod)] <- NA
                                     return(rank(-sim_row) - 1)}, # max similarity will be ranked 1st.
                                mc.cores = threads))

    }
    #self_rank = mapply(function(row, i) row[i], as.data.frame(rank_matrix), reflection_idx)
    self_rank = unlist(mclapply(1:nrow(rank_matrix), function(i) rank_matrix[i, reflection_idx[i]], 
                                mc.cores = threads))
    return(self_rank)
}

calculate_self_dist_rank <- function(embedding, cell_identity, modality, threads = 20, 
                                    dissim = c("Euclidean", "dot product")){
    require(foreach)
    R = 20
    cl <- parallel::makeForkCluster(threads)
    doParallel::registerDoParallel(cl)
    # Samples 5% of the cell identity to check the neighborhood ranks.
    sample_size = min(2000, length(cell_identity) %/% R)
    self_ranks.l = foreach(r=1:R, .combine = 'c') %dopar% {
        list(sampled_self_rank(embedding, cell_identity, modality, sample_size, seed = r) / sample_size)
    }
    
    parallel::stopCluster(cl)
    return(sapply(self_ranks.l, median))
    #return(sapply(self_ranks.l, median))
}

calculate_self_dist_rank_all <- function(embedding, cell_identity, modality, threads = 20,
                                        dissim = c("Euclidean", "dot product")){
    sample_size = length(cell_identity) %/% 2
    self_ranks.l = sampled_self_rank(embedding, cell_identity, modality, sample_size, seed = 1, threads = threads, dissim = dissim) / sample_size

    return(self_ranks.l)
                #return(sapply(self_ranks.l, median))
}


calculate_anchoring_dist <- function(embedding, cell_identity, k = 200, summary_stat = median, 
                                    dissim = c("Euclidean", "dot product")){
    reflection_idx = get_reflection_idx(cell_identity)
    my_reflection = embedding[reflection_idx,]

    if (dissim == "Euclidean"){
        dist_self = calculate_pdist(embedding, my_reflection)

        # Get normalization factor: mean distance to kNN (k = 50)
        #require(FNN)
        #res = get.knn(embedding, k)
        #norm_factor = apply(res$nn.dist, 1, mean)

        # Normalize by randomly sampled distance
        nsample_pairs = (length(cell_identity)/2)%/%10
        message(paste0("Normalizing for mean distance of ", nsample_pairs, " random pairs"))
        dists = sapply(1:nsample_pairs, 
                        function(i) {
                            pair = sample(cell_identity, 2); 
                            calculate_pdist(embedding[pair[1],], embedding[pair[2],])
                        }
                      )
        norm_factor = mean(dists)

        anchoring_dist = dist_self / norm_factor 
    
    } else if (dissim == "dot product") {
        dist_matrix = get_dotsim_dist_matrix(embedding)                            
        dist_self = sapply(1:nrow(dist_matrix), function(i) dist_matrix[i, reflection_idx[i]])
        
        # Randomly sample similarities
        nsample_pairs = (length(cell_identity)/2)%/%10
        message(paste0("Normalizing for mean similarity of ", nsample_pairs, " random pairs"))
        dists = sapply(1:nsample_pairs,
                        function(i) {
                            set.seed(i)
                            pair = sample(cell_identity, 2);
                            dist_matrix[pair[1], pair[2]]
                        }
                      )
        norm_factor = mean(dists)

        anchoring_dist = dist_self / norm_factor
 
    }
    
    
    return(summary_stat(anchoring_dist))
}

get_n_clusters <- function(embedding, nClusters, max_steps = 20, dissim){
    require(Seurat)
    if (dissim == "Euclidean"){
        nn_graph = FindNeighbors(embedding)
    } else {
        dist_matrix = get_dotsim_dist_matrix(embedding)
        nn_graph = FindNeighbors(embedding, distance.matrix = )
    }
    steps = 0
    current_min = 0
    current_max = 3
    while(steps < max_steps){
        print(paste0('step ', steps))
        current_resolution = current_min + ((current_max - current_min)/2)
        clusters = FindClusters(nn_graph$snn, resolution = current_resolution)[,1]
        current_n = length(unique(clusters))
        print(paste0("got ", current_n, " at resolution ", current_resolution))
        if (current_n > nClusters) current_max = current_resolution
        else if (current_n < nClusters) current_min = current_resolution
        else{return(clusters)}
        steps = steps + 1
    }
    print("Cannot find the number of clusters")
    print(paste0("Clustering solution from last iteration is used:", current_n, " at resolution ", current_resolution))
    return(clusters)
}

calculate_ARI <- function(embedding, cell_identity, nClusters = 20, 
                          dissim = c("Euclidean", "dot product")){
    require(Seurat)
    clusters = get_n_clusters(embedding, nClusters = 20, dissim = dissim) #res = 0.8 by default
    in_same_celltype = sapply(levels(cell_identity), function(lvl){
                          cell_types = clusters[which(cell_identity == lvl)]
                          if( length(cell_types) != 2) {print(lvl); stop()}
                          return( as.numeric(cell_types[1] == cell_types[2]) )})
    print(str(in_same_celltype))
    return(in_same_celltype)
}

get_summarisable_metrics <- function(embedding, cell_identity, summary_stat, dissim){
    summary_sil = calculate_silhouette(embedding, embedding$cell.ft, summary_stat = summary_stat, 
                                        dissim = dissim)
    summary_dist_rank = calculate_self_dist_rank(embedding, embedding$cell.ft, summary_stat = summary_stat, dissim = dissim)
    summary_anchor_dist = calculate_anchoring_dist(embedding, embedding$cell.ft, summary_stat = summary_stat, dissim = dissim)
    res.df = data.frame(Silhouette = summary_sil, 
                        Anchoring_dist_rank = summary_dist_rank, 
                        Anchoring_dist = summary_anchor_dist, 
                        summary_stat = as.character(substitute(summary_stat)))
    message("Done calculating summarisable metrics")
    print(res.df)
    return(res.df)  
}

get_value_metrics <- function(embedding, cell_identity, summary_stat, dissim){
    metrics = get_summarisable_metrics(embedding, cell_identity, summary_stat, dissim = dissim)
    metrics$ARI = calculate_ARI(embedding, dissim = dissim)
    metrics$summary_method = as.character(substitute(summary_stat))
    return(metrics)
} 

write_value_metrics <- function(embedding, cell_identity, dissim){
    medians = get_value_metrics(embedding, cell_identity, median, dissim = dissim)
    means = get_value_metrics(embedding, cell_identity, mean, dissim = dissim)
    write.table(rbind(medians, means), file = paste0(output_path,"/summary_metrics.txt"), 
                quote = F, row.names = T, col.names= T, sep = "\t")
    message("wrote summary metrics")
}

write_cell_metrics <- function(embedding, cell_identity, dissim){
    res.df = get_summarisable_metrics(embedding, cell_identity, identity, dissim = dissim)
    write.table(res.df, file = paste0(output_path,"/per_cell_metrics.txt"), 
                quote = F, row.names = T, col.names = T, sep = "\t")
    message("wrote per cell metrics")
}

