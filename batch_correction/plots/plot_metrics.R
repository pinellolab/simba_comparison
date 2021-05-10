args = commandArgs(trailingOnly = TRUE)

metrics_file_path = args[1]
output_plot_path = args[2]

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

library(cowplot)
theme_set(theme_cowplot(font_size=20))

cbPalette <- c("#000000", "#999999", "#D55E00", "#009E73",
                               "#E69F00", "#0072B2") #"#009E73","#F0E442", "#CC79A7", "#56B4E9", 
metrics <- read.table(metrics_file_path, head = T)
metrics$batch_sub = 1 - metrics$batch
metrics$celltype_sub = 1 - metrics$celltype

metrics$method = factor(metrics$method, levels = c("Raw", "Raw-pp", "Seurat3", "LIGER", "Harmony", "SIMBA"))

asw = subset(metrics, metric %in% c("ASW"))
ari = subset(metrics, metric %in% c("ARI"))
lisi = subset(metrics, metric %in% c("LISI_40"))

get_limit <- function(v){
    range_len = max(v) - min(v)
    return(c(min(v) - 0.1*range_len, max(v) + 0.1*range_len))
}

p1 <- ggplot(asw, aes(celltype, batch_sub, col = method, label = method)) + geom_point(size = 1) + theme_classic(base_size = 15) + xlab("ASW cell type") + ylab("1 - ASW batch") + scale_color_manual(values =  cbPalette) + theme(legend.position = "none") + geom_text_repel(force = 50, size = 5) + scale_x_continuous(limits = get_limit(asw$celltype)) + scale_y_continuous(limits = get_limit(asw$batch_sub))
p2 <- ggplot(ari, aes(celltype, batch_sub, col = method, label = method)) + geom_point(size = 1) + theme_classic(base_size = 15) + xlab("ARI cell type") + ylab("1 - ARI batch") + scale_color_manual(values =  cbPalette) + theme(legend.position = "none") + geom_text_repel(force = 50, force_pull = 0.01, size = 5) + scale_x_continuous(limits = get_limit(ari$celltype)) + scale_y_continuous(limits = get_limit(ari$batch_sub))
p3 <- ggplot(lisi, aes(celltype_sub, batch, col = method, label = method)) + geom_point(size = 1) + theme_classic(base_size = 15) + xlab("1 - cLISI (cell type)") + ylab("iLISI batch") + scale_color_manual(values =  cbPalette) + theme(legend.position = "none") + geom_text_repel(size = 5) + scale_x_continuous(limits = get_limit(lisi$celltype_sub)) + scale_y_continuous(limits = get_limit(lisi$batch))

pdf(output_plot_path, height = 3, width = 9, useDingbats = FALSE)
plot_grid(p1, p2, p3, ncol = 3)
dev.off()
