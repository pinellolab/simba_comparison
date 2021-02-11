args = commandArgs(trailingOnly = TRUE)
output_dir = args[length(args)]
dataset = args[1]
input_files = args[2:(length(args)-1)]

split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))

input_files.l = lapply(input_files, function(path){
                       df = read.table(path, head = T, stringsAsFactors = F)
                       s = split_path(path)
                       software = s[3]
                       dataset = s[2]
                       df$software = software
                       df$dataset = dataset
                        return(df)
})

combined.df = do.call(rbind, input_files.l)
write.table(combined.df, paste0(output_dir, "/metrics_collected_", dataset, ".txt"), quote = F, sep = "\t", col.names = T, row.names = T)

input.files.cell.l = lapply(input_files, function(path){
                       s = split_path(path)
                       software = s[3]
                       dataset = s[2]
                       df = read.table(paste0(software, "/", dataset, "/per_cell_metrics.txt"), sep = "\t", head = T)
                       df$software = software
                       df$dataset = dataset
                       return(df)

})
combined.cell.df = do.call(rbind, input.files.cell.l)
write.table(combined.df, paste0(output_dir, "/per_cell_metrics_", dataset, ".txt"), quote = F, sep = "\t", col.names = T, row.names = T)
