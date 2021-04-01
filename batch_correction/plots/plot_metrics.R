args = commandArgs(trailingOnly = TRUE)

metrics_file_path = args[1]
output_plot_path = args[2]

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

cbPalette <- c("#000000", "#999999", "#D55E00", "#009E73",
                               "#E69F00", "#0072B2") #"#009E73","#F0E442", "#CC79A7", "#56B4E9", 
metrics <- read.table(metrics_file_path, head = T)
metrics$batch_sub = 1 - metrics$batch
metrics$celltype_sub = 1 - metrics$celltype

metrics$method = factor(metrics$method, levels = c("Raw", "Raw-pp", "Seurat3", "LIGER", "Harmony", "SIMBA"))

p1 <- ggplot(subset(metrics, metric %in% c("ASW")), aes(celltype, batch_sub, col = method, label = method)) + geom_point(size = 3) + theme_classic(base_size = 20) + xlab("ASW cell type") + ylab("1 - ASW batch") + scale_color_manual(values =  cbPalette) + theme(legend.position = "none") + geom_text_repel(force = 50)
p2 <- ggplot(subset(metrics, metric %in% c("ARI")), aes(celltype, batch_sub, col = method, label = method)) + geom_point(size = 3) + theme_classic(base_size = 20) + xlab("ARI cell type") + ylab("1 - ARI batch") + scale_color_manual(values =  cbPalette) + theme(legend.position = "none") + geom_text_repel(force = 50, force_pull = 0.01)
p3 <- ggplot(subset(metrics, metric %in% c("LISI_40")), aes(celltype_sub, batch, col = method, label = method)) + geom_point(size = 3) + theme_classic(base_size = 20) + xlab("1 - cLISI (cell type)") + ylab("iLISI batch") + scale_color_manual(values =  cbPalette) + theme(legend.position = "none") + geom_text_repel()

pdf(output_plot_path, height = 5, width = 15)
cowplot::plot_grid(p1, p2, p3, ncol = 3)
dev.off()
