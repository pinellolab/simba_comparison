args = commandArgs(trailingOnly = T)

embedding_file = args[1]
output_file = args[2]
if (length(args) > 2) {nDims = as.numeric(args[3])} else {nDims = 20}
tbl = read.table(embedding_file, sep = "\t", row.names = 1, head = T)
embedding = tbl[,1:(ncol(tbl)-2)]

emb.pca = prcomp(as.matrix(embedding), center = T, scale. = F, rank. = nDims)
out.df = as.data.frame(emb.pca$x)

cid = which(colnames(tbl) %in% c("celltype", "cell_type", "CellType"))
bid = which(colnames(tbl) %in% c("batch", "batchlb"))

if(length(cid) !=1){print(cid); stop()}
stopifnot(length(bid) == 1)

out.df$cell_type = tbl[,cid]
out.df$batch = tbl[,bid]

write.table(file = output_file, out.df, quote = F, sep = "\t", row.names= T, col.names = T)
