args = commandArgs(trailingOnly = TRUE)

umap_file = args[1]
output_file = args[2]
palette = args[3]
if (length(args) < 4){ labels = FALSE } else{
    if (args[4] == "TRUE") {labels = TRUE
    } else {
        labels = FALSE}
}


tbl <- read.table(umap_file, sep = "\t", head = T, stringsAsFactors = FALSE)
tbl$cell_type[tbl$cell_type == "MHC_class_II"] <- 'MHC class II'
require(ggplot2)

tbl$cell_type = factor(tbl$cell_type)
tbl$batch = factor(tbl$batch)
p1 = Seurat:::SingleDimPlot(tbl, dims = c("UMAP_1", "UMAP_2"), col.by = "batch", raster = TRUE)
p2 = Seurat:::SingleDimPlot(tbl, dims = c("UMAP_1", "UMAP_2"), col.by = "cell_type", raster = TRUE)

if (palette == "pancreas"){
    #palette_celltype={'TAC-1':'#F8D856', 'TAC-2':'#F1B044', 'IRS':'#779ba1',               'Medulla':'#897a74','Hair Shaft-cuticle.cortex':"#d6a780"}
    pal = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', "#aa40fc",  '#8c564b', '#e377c2', '#b5bd61', '#17becf', '#aec7e8', 
            '#ffbb78', '#98df8a', '#ff9896',  '#c5b0d5', '#c49c94')
    tbl$cell_type = factor(tbl$cell_type, levels = c('MHC class II', 'acinar', 'alpha', 'beta', 'delta', 'ductal', 'endothelial', 'epsilon', 'gamma', 'macrophage', 'mast',
                                                     'mesenchymal', 'schwann', 'stellate', 't_cell' ))
    pal_batch = c('#4c72b0', '#dd8452', '#55a868', '#c44e52', '#8172b3')
    p1 = Seurat:::SingleDimPlot(tbl, dims = c("UMAP_1", "UMAP_2"), col.by = "batch", raster = TRUE)
    p2 = Seurat:::SingleDimPlot(tbl, dims = c("UMAP_1", "UMAP_2"), col.by = "cell_type", raster = TRUE)

    p1 = p1 + scale_color_manual(values = pal_batch)
    p2 = p2 + scale_color_manual(values = pal)

} else if (palette == "murine-atlas") {
    pal = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78')
    tbl$cell_type = factor(tbl$cell_type, levels = c("B-cell", "Dendritic", 'Endothelial', 'Epithelial', 'Macrophage', 'Monocyte', 'NK', 'Neutrophil', 'Smooth-muscle', 'Stromal', 'T-cell'))
    pal_batch = c('#4c72b0', '#dd8452')
    p1 = p1 + scale_color_manual(values = pal_batch)
    p2 = p2 + scale_color_manual(values = pal)
}

if (! labels) {
    p1 = p1 + theme(legend.position = "none")
    p2 = p2 + theme(legend.position = "none")
}

pdf(output_file, width = 10, height = 5)

cowplot::plot_grid(p1, p2)
dev.off()


