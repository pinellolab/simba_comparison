library(dplyr)

args = commandArgs(trailingOnly = TRUE)
nArgs = length(args)
method_use = args[1]
ARI_output_path = args[2]
ASW_output_path = args[3]
LISI_output_paths = args[4:(nArgs-3)]
output_dir = args[nArgs-2]
emb_type = args[nArgs-1]
dissim_metric = args[nArgs]

if ( ! grepl("ARI", ARI_output_path) ) stop(paste0("First input not ARI:", ARI_output_path))
if ( ! grepl("ASW", ASW_output_path) ) stop("Second input not ASW")

ARI_output <- read.table(ARI_output_path, head = T, stringsAsFactors = F)
ARI_median <- ARI_output[grepl("median", ARI_output$use_case), c("use_case", "ari_batch", "ari_celltype")]
stopifnot(length(ARI_median$use_case) == 1)
ARI_median[1, "use_case"] = "ARI"
colnames(ARI_median) <- c("metric", "batch", "celltype")

ASW_output <- read.table(ASW_output_path, head = T, row.names = 1, stringsAsFactors = F)
ASW_median <- ASW_output[grepl("median", ASW_output$method_use), c("method_use", "asw_batch", "asw_celltype")]
stopifnot(length(ASW_median$method_use) == 1)
ASW_median[1, "method_use"] = "ASW"
colnames(ASW_median) <- c("metric", "batch", "celltype")

collected_metrics = bind_rows(list(ARI_median, ASW_median))
for (LISI_output_path in LISI_output_paths){
    LISI_output <- read.table(LISI_output_path, head = T, stringsAsFactors = F)
    LISI_dat <- LISI_output[LISI_output$methods_use == method_use, c("methods_use", "iLISI_median", "cLISI_median")]
    LISI_dat[1, "methods_use"] = paste0("LISI_", LISI_output$plx[1])
    colnames(LISI_dat) <- c("metric", "batch", "celltype")
    collected_metrics = bind_rows(list(collected_metrics, LISI_dat))
}

write.table(collected_metrics, file = paste0(output_dir, "/", method_use, "_metrics_collected_", emb_type, "_", dissim_metric, ".txt"), row.names = F, quote = F)
