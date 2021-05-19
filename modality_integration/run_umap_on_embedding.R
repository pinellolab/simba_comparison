args = commandArgs(trailingOnly = T)

embedding_file = args[1]
output_file = args[2]

tbl = read.table(embedding_file, sep = "\t", head = T, row.names= 1)
embedding = tbl[,1:(ncol(tbl)-2)]

emb.umap = uwot::umap(embedding, n_components = 2, 
                     n_neighbors = 30, 
                     min_dist = 0.3, 
                     metric = "cosine")
out.df = as.data.frame(emb.umap)
rownames(out.df) = rownames(embedding)
colnames(out.df) = c("UMAP_1", "UMAP_2")

cid = which(colnames(tbl) %in% c("celltype", "cell_type", "CellType"))
bid = which(colnames(tbl) %in% c("batch", "batchlb"))

if(length(cid) !=1){print(cid); stop()}
stopifnot(length(bid) == 1)

out.df$cell_type = tbl[,cid]
out.df$batch = tbl[,bid]
out.df$modality = sapply(rownames(tbl), function(s) strsplit(s, split = ".", fixed = TRUE)[[1]][1])
write.table(file = output_file, out.df, quote = F, sep = "\t", row.names= T, col.names = T)
