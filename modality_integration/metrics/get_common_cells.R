args = commandArgs(trailingOnly = T)

input_files = args[1:(length(args) - 1)]
output_file = args[length(args)]

cells.l = lapply(input_files, function(path){ rownames(read.table(path)) })
common_cells = Reduce(intersect, cells.l)

message(paste0(length(common_cells), " common cells saved in ", output_file))
head(common_cells)
write.table(file = output_file, common_cells, col.names = F, quote = F, row.names = F)
