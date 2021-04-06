source("metrics/metric_utils.R")

args = commandArgs(trailingOnly = TRUE)

metric = args[1]
embedding_file = args[2]
cells_file = args[3]
software = args[4]
output_path = args[5]
nThreads = as.numeric(args[6])

stopifnot(metric %in% c("Silhouette", "Anchoring_dist_rank", "Anchoring_dist", "Cluster_agreement", "Anchoring_dist_rank_full"))

cells_use = read.table(cells_file, head = F, stringsAsFactors = F)$V1
head(cells_use)
embedding = read.table(embedding_file)
embedding = embedding[cells_use, ]
nDims = ncol(embedding)-2
cell_identity = factor(embedding$batch)
embedding = embedding[,1:nDims]
nCells = length(cell_identity)/2

if (software == "SIMBA-dot") {
    dissim = "dot product"
} else {
    dissim = "Euclidean"
}

# get modality
modality = sapply(strsplit(rownames(embedding), split = ".", fixed = TRUE), '[', 1)

if (metric == "Silhouette") {
    res = calculate_silhouette(embedding, cell_identity, summary_stat = identity, dissim = dissim)
} else if (metric == "Anchoring_dist_rank") {
   res = calculate_self_dist_rank(embedding, cell_identity, modality = modality, threads = 1, dissim = dissim)
} else if (metric == "Anchoring_dist") {
    res = calculate_anchoring_dist(embedding, cell_identity, summary_stat = identity, dissim = dissim)
}  else if (metric == "Cluster_agreement") {
    res = calculate_ARI(embedding, cell_identity, dissim = dissim)
}  else if (metric == "Anchoring_dist_rank_full") {
   res = calculate_self_dist_rank_all(embedding, cell_identity, modality = modality, threads = nThreads, dissim = dissim)
} 

res.df = data.frame(metric = metric, value = res, mean = mean(res), median = median(res))
rownames(res.df) = names(res)
write.table(res.df, paste0(output_path, "/metric_", metric, ".txt"), sep = "\t", quote = F, col.names = F, row.names= T)
message(paste0("Done producing ", metric))
