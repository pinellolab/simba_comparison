library(GenomicRanges)
library(rtracklayer)

args = commandArgs(trailingOnly = TRUE)
if(length(args) != 3){
    message("Usage: Rscript make_promoter_regions.R gene_file.bed promoter_file.bed output_file_name")
}
gene_file = args[1]
upstream_length = as.numeric(args[2])
output_file = args[3]

genes <- import(gene_file)
pms = restrict(promoters(genes, upstream = upstream_length, downstream = 0), start = 1, keep.all.ranges = TRUE)

df.pm <- data.frame(seqnames = seqnames(pms), starts = start(pms) - 1, ends = end(pms), names = pms$name, scores = c(rep(".", length(pms))), strans = strand(pms))
write.table(df.pm, file = output_file, quote = F, sep = "\t", row.names = F, col.names = F)

