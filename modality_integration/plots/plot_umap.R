args = commandArgs(trailingOnly = TRUE)

umap_file = args[1]
output_file = args[2]
palette = args[3]
if (length(args) < 4){ labels = FALSE } else{
    if (args[4] == "TRUE") {labels = TRUE
    } else {
        labels = FALSE}
}


tbl <- read.table(umap_file, sep = "\t", head = T)
require(ggplot2)

pal_batch = c('#4c72b0', '#dd8452')

tbl$cell_type = factor(tbl$cell_type)
if (length(tbl$modality) == 0){tbl$modality = sapply(rownames(tbl), function(s) strsplit(s, split = ".", fixed = TRUE)[[1]][1])}
tbl$modality = factor(tbl$modality)

print(str(tbl$modality))
tbl <- tbl[sample(1:nrow(tbl)),]

p1 = Seurat:::SingleDimPlot(tbl, dims = c("UMAP_1", "UMAP_2"), col.by = "modality", raster = TRUE, order = ) + scale_color_manual(values = pal_batch)
p2 = Seurat:::SingleDimPlot(tbl, dims = c("UMAP_1", "UMAP_2"), col.by = "cell_type", raster = TRUE, order = )

if (palette == "mouse-skin" || palette == "mouse-skin-subset"){
    #palette_celltype={'TAC-1':'#F8D856', 'TAC-2':'#F1B044', 'IRS':'#779ba1',               'Medulla':'#897a74','Hair Shaft-cuticle.cortex':"#d6a780"}
    pal = c('#F8D856', '#F1B044', '#C37777', '#897a74', "#d6a780")

    tbl$cell_type = factor(tbl$cell_type, levels = c("TAC-1", "TAC-2", "IRS", "Medulla", "Hair_Shaft-cuticle.cortex"))
    p2 = Seurat:::SingleDimPlot(tbl, dims = c("UMAP_1", "UMAP_2"), col.by = "cell_type", raster = TRUE)
    p2 = p2 + scale_color_manual(values = pal)

} else if (palette == "pbmc") {
    pal = c('#4c72b0', '#dd8452', '#55a868', '#c44e52', '#8172b3', '#937860', '#da8bc3', '#8c8c8c', '#ccb974', '#64b5cd', '#a1c9f4', '#ffb482', '#8de5a1', '#ff9f9b', '#d0bbff', '#debb9b', '#fab0e4', '#cfcfcf', '#fffea3')
    p2 = p2 + scale_color_manual(values = pal)
}

if (! labels) {
    p1 = p1 + theme(legend.position = "none")
    p2 = p2 + theme(legend.position = "none")
}

pdf(output_file, width = 10, height = 5, useDingbats = FALSE)

cowplot::plot_grid(p1, p2)
dev.off()


