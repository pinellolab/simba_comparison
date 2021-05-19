library(plyr)
library(dplyr)
library(ggplot2)
library(egg)

args = commandArgs(trailingOnly = TRUE)
metric_file = args[1]
output_pdf = args[2]

methods_use = c("Raw", "Seurat3", "LIGER", "SIMBA")
cbPalette <- c("#999999",  "#D55E00","#009E73", "#0072B2","#E69F00", "#56B4E9", "#F0E442", "#CC79A7")
res = read.table(metric_file)

colnames(res)[6] <- "method"
res = res %>% subset(method %in% methods_use)
res$method <- factor(res$method, levels = methods_use)

median_cl_boot <- function(x, conf = 0.95) {
    lconf <- (1 - conf)/2
    uconf <- 1 - lconf
    require(boot)
    bmedian <- function(x, ind) median(x[ind])
    bt <- boot(x, bmedian, 3000)
    bb <- boot.ci(bt, type = "perc")
    data.frame(y = median(x), ymin = quantile(bt$t, lconf), ymax = quantile(bt$t, 
        uconf))
}

require(boot)
boot_CIs <- function(df, conf = 0.95){
    lconf <- (1 - conf)/2
    uconf <- 1 - lconf
    bmean = function(x, ind) mean(x[ind,]$value)
    bootobject = boot(df, bmean, R = 3000)
    boot_ci = boot.ci(bootobject, type = "perc", conf = conf)
    CI_df = data.frame(y = mean(df$value), ymin = quantile(bootobject$t, lconf), ymax = quantile(bootobject$t, 
            uconf))
    return(CI_df)
}

boot_df.l = lapply(levels(res$method), function(m) {
    df = subset(res, metric == "Cluster_agreement" & method == m)
    return(boot_CIs(df))
})


names(boot_df.l) <- levels(res$method)
boot_df = do.call(rbind, boot_df.l)
boot_df$method = factor(rownames(boot_df), levels = methods_use)

options(repr.plot.width = 20, repr.plot.height = 5)
df = subset(res, metric == "Anchoring_dist_rank_full") %>% arrange(value) %>% group_by(method) %>% mutate(rank = row_number())
p1 = ggplot(df, aes(x = rank, y = value, color = method)) + ylab("Anchoring distance rank") + geom_line() + theme_classic(base_size = 18) + scale_color_manual(values = cbPalette) + xlab("cells")
p2 = ggplot(subset(res, metric == "Anchoring_dist"), aes(x = method, y = value, fill = method)) + ylab("Anchoring distance") + geom_violin() + geom_boxplot(width = 0.1) + theme_classic(base_size = 18) + scale_y_log10() + scale_fill_manual(values = cbPalette) + theme(legend.position = "none") + stat_summary(fun.data = median_cl_boot, geom = "errorbar", colour = "red") + stat_summary(fun = median, geom = "point", colour = "red")
p3 = ggplot(subset(res, metric == c("Silhouette")), aes(x = method, y = value, fill = method)) + ylab("Silhouette index") + geom_violin() + geom_boxplot(width = 0.1, alpha = 0, position = 'dodge') + theme_classic(base_size = 18) + scale_fill_manual(values = cbPalette) + theme(legend.position = "none")+ stat_summary(fun.data = median_cl_boot, geom = "errorbar", colour = "red")
p4 = ggplot(boot_df, aes(x = method, y = y, color = method)) + ylab("Fraction in same cluster") + geom_point(shape = 4, size = 5, stroke = 2) + theme_classic(base_size = 18) + scale_color_manual(values = cbPalette) + theme(legend.position = "none") + geom_errorbar(aes(ymin=ymin, ymax=ymax))

pdf(output_pdf, width = 20, height = 5, useDingbats = FALSE)
ggarrange(p1, p2, p3, p4, ncol = 4)
dev.off()

