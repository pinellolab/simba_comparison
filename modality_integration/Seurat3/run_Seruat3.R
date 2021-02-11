# https://satijalab.org/seurat/archive/v3.0/atacseq_integration_vignette.html

library(Seurat)
library(ggplot2)
library(patchwork)
library(Signac)

args = commandArgs(trailingOnly = TRUE)
annotation_file= "modality_integration/data/Homo_sapiens.GRCh37.82.gtf.gz"

main <- function(rds_path){
data <- readRDS(args[1])
dataset_name = args[2]
peaks <- data[['Peaks']]
genes = data[['Gene Expression']]
    atac.sobj = process_atac(peaks)
    rna.sobj = process_rna(genes)
}

process_atac <- function(peak_matrix, dataset_name, 
                        nCount_ATAC_thres = 5000,
                        nSV = 50){
# create a gene activity matrix from the peak matrix and GTF, using chromosomes 1:22, X, and Y.
# Peaks that fall within gene bodies, or 2kb upstream of a gene, are considered
    activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peak_matrix, 
                                                annotation.file = annotation_file, 
                                                seq.levels = c(1:22, "X", "Y"), upstream = 2000, verbose = TRUE)

    pbmc.atac <- CreateSeuratObject(counts = peak_matrix, assay = "ATAC", project = paste0(dataset_name, "_ATAC"))
    pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)

    pbmc.atac <- subset(pbmc.atac, subset = nCount_ATAC > nCount_ATAC_thres)
    pbmc.atac$tech <- "atac"

    DefaultAssay(pbmc.atac) <- "ACTIVITY"
    pbmc.atac <- FindVariableFeatures(pbmc.atac)
    pbmc.atac <- NormalizeData(pbmc.atac)
    pbmc.atac <- ScaleData(pbmc.atac)

    DefaultAssay(pbmc.atac) <- "ATAC"
    VariableFeatures(pbmc.atac) <- names(which(Matrix::rowSums(pbmc.atac) > 100))
    pbmc.atac <- RunLSI(pbmc.atac, n = nSV, scale.max = NULL)
    pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 1:nSV)
    return(pbmc.atac)

}

process_rna <- function(gene_matrix, dataset_name, 
                        min.cells = 3, 
                        min.features = 200, 
                        percent.mt.thres = 20, 
                        nCount_RNA_max = 25000,
                        nHVGs = 2000, 
                        nCount_RNA_min = 1000, 
                        



















































