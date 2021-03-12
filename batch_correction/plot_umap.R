args = commandArgs(trailingOnly = TRUE)

umap_file = args[1]
output_file = args[2]

tbl <- read.table(umap_file, sep = "\t")
require(ggplot2)
pdf(output_file, width = 10, height = 5)
cowplot::plot_grid(Seurat:::SingleDimPlot(tbl, dims = c("UMAP_1", "UMAP_2"), col.by = "batch") + theme(legend.position = "none"), 
                   Seurat:::SingleDimPlot(tbl, dims = c("UMAP_1", "UMAP_2"), col.by = "cell_type") + theme(legend.position = "none"))
dev.off()
