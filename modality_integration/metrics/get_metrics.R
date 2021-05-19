source("metrics/metric_utils.R")

args = commandArgs(trailingOnly = TRUE)

embedding_file = args[1]
nDims = as.numeric(args[2])
output_path = args[3]
nThreads = as.numeric(args[4])
load_files = TRUE
save_files = TRUE



embedding = read.table(embedding_file)
embedding$cell.ft = factor(unlist(mclapply(strsplit(rownames(embedding), split = ".", fixed = TRUE), function(cid) paste0(cid[2:length(cid)], collapse = "."), mc.cores = 20)))

if (load_files & file.exists(paste0(output_path, "/.dist.m.rds")) & file.exists(paste0(output_path, "/.dist.rank.m.rds"))){
    dist.m = readRDS(paste0(output_path, "/.dist.m.rds"))
    dist.rank.m = readRDS(paste0(output_path, "/.dist.rank.m.rds")) 
} else {
    message("Producing distance matrices")
    dist.m = as.matrix(parDist(as.matrix(embedding[,1:nDims]), method = "euclidean", threads = nThreads))
    dist.rank.m = apply(dist.m, 1, FUN=rank) - 1
    if (save_files){
    saveRDS(dist.m, paste0(output_path, "/.dist.m.rds"))
    saveRDS(dist.rank.m, paste0(output_path, "/.dist.rank.m.rds"))}
}

message("Done producing distance matrices")


write_value_metrics(dist.m, dist.rank.m, embedding, embedding$cell.ft)
write_cell_metrics(dist.m, dist.rank.m, embedding, embedding$cell.ft)
