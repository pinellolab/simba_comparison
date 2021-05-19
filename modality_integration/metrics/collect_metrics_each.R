args = commandArgs(trailingOnly = TRUE)
output_dir = args[length(args)]
input_files = args[1:(length(args)-1)]

split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))

input_files.l = lapply(input_files, function(path){
                       df = read.table(path, head = F, stringsAsFactors = F)
                       colnames(df) = c("cell", "metric", "value", "mean", "median") 
                        return(df)
})

combined.df = do.call(rbind, input_files.l)
write.table(combined.df, paste0(output_dir, "/metrics_collected.txt"), quote = F, sep = "\t", col.names = T, row.names = F)
message("Metrics collected in ", output_dir)
