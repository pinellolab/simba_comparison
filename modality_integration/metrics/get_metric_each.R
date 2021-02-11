source("metrics/metric_utils.R")

args = commandArgs(trailingOnly = TRUE)

metric = args[1]
embedding_file = args[2]
nDims = as.numeric(args[3])
output_path = args[4]
load_files = TRUE
save_files = TRUE

stopifnot(metric %in% c("Silhouette", "Anchoring_dist_rank", "Anchoring_dist", "Graph_connectivity", "Celltype_agreement"))


embedding = read.table(embedding_file)
embedding$cell.ft = factor(unlist(lapply(strsplit(rownames(embedding), split = ".", fixed = TRUE), function(cid) paste0(cid[2:length(cid)], collapse = "."))))

load_or_create <- function(tmp_rds){
    if (load_files & file.exists(paste0(output_path, tmp_rds))) dist.m = readRDS(paste0(output_path, tmp_rds))
    else {
        dist.m = as.matrix(dist(embedding[,1:nDims], method = "euclidean"))
        if (save_files){
                saveRDS(dist.m, paste0(output_path, tmp_rds))
        }
    }

    message("Done producing distance matrices")
    return(dist.m)
}




if (metric == "Silhouette") {
    dist.m = load_or_create("/.dist.m.rds")
    res = calculate_silhouette(dist.m, embedding$cell.ft, summary_stat = identity)
} else if (metric == "Anchoring_dist_rank") {
   dist.rank.m = load_or_create("/.dist.rank.m.rds")
   res = calculate_self_dist_rank(dist.rank.m, embedding$cell.ft, summary_stat = identity)
} else if (metric == "Anchoring_dist") {
    dist.rank.m = load_or_create("/.dist.rank.m.rds")
    res = calculate_self_dist_rank(dist.rank.m, embedding$cell.ft, summary_stat = identity)
} else if (metric == "Graph_connectivity") {
    dist.rank.m = load_or_create("/.dist.rank.m.rds")
    res = calculate_graph_connectivity(dist.rank.m, embedding$cell.ft, k = 50)
} else if (metric == "Celltype_agreement") {
    res = calculate_ARI(embedding)
}

res.df = data.frame(metric = metric, value = res, mean = mean(res), median = median(res))
rownames(res.df) = names(res)
write.table(res.df, paste0(output_path, "/metric_", metric, ".txt"), sep = "\t", quote = F, col.names = F, row.names= T)
message(paste0("Done producing ", metric))
