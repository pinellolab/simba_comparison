library(dplyr)
args = commandArgs(trailingOnly = TRUE)

out_dir = args[length(args)-1]
nPCs = args[length(args)]
table_paths = args[1:(length(args)-2)]
tables <- lapply(table_paths, function(table_path) read.table(table_path, head = T))
print(tables)

get_method_name <- function(file_path){
    return(strsplit(basename(file_path), split = "_metrics_collected")[[1]][1])
}

method_names <- sapply(table_paths, get_method_name)
names(tables) = method_names

cat_table = bind_rows(tables, .id = "method")
write.table(cat_table, file = paste0(out_dir, "/metrics_collected_PC", nPCs, ".txt"), quote = F, sep = "\t", row.names = F)
                    
