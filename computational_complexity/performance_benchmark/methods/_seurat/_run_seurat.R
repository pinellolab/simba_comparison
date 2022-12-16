

# -- parse incoming args: ----------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)

h5ad_path = as.character(args[1])
obs_path = as.character(args[2])
var_path = as.character(args[3])
write_dir = as.character(args[4])
project_name = as.character(args[5])
N = as.integer(args[6])

# -- load libraries: ---------------------------------------------------------------------

library(dplyr)
library(Seurat)
library(stringr)
library(patchwork)
library(SeuratData)
library(SeuratDisk)

# -- function to collect memory usage: ---------------------------------------------------

track_mem <- function(mem_info_start, fname) {

    mem_info_end = gc()

    mem_df <- as.data.frame(c(
    "mem_start_current"=mem_info_start[9],
    "mem_start_max"=mem_info_start[10],
    "mem_end_current"=mem_info_end[9],
    "mem_end_max"=mem_info_end[10]
    ))

    write.csv(mem_df, paste(fname, ".csv", sep=""))
}

time_dict <- c()

# -- convert .h5ad -> .h5Seurat: ---------------------------------------------------------

h5seurat_path = paste(str_split(h5ad_path, ".h5ad")[[1]][1], ".h5seurat", sep="")

print(h5seurat_path)

if (file.exists(h5seurat_path)) {
    # do nothing
    } else {
    Convert(h5ad_path, dest = "h5seurat", overwrite = TRUE)
    }


# -- read h5 and .csv from adata to Seurat: ----------------------------------------------

dataset <- LoadH5Seurat(h5seurat_path)
obs_df = read.csv(obs_path)
var_df = read.csv(var_path)


# -- annotate adata with .obs and .var: --------------------------------------------------

obs_cols = colnames(obs_df)
var_cols = colnames(var_df)

for (i in obs_cols){
    dataset[[i]] = obs_df[[i]]
}

for (i in var_cols){
    dataset[["RNA"]][[i]] = var_df[[i]]
}

# -- downsample: -------------------------------------------------------------------------
if (N > 0) {
    dataset = subset(x = dataset, downsample = N)
} else {
}

counts = GetAssayData(object = dataset, slot="counts")
colnames(counts) = dataset$barcodes
row.names(counts) = dataset[['RNA']][[]]$gene_names


# -- begin timing and memory observation: ------------------------------------------------
mem_info_pipeline_init = gc(reset=TRUE)
mem_info_start = gc(reset=TRUE)
time_dict['t0'] = t0 = Sys.time()
# ----------------------------------------------------------------------------------------
print(time_dict)

t_ref = Sys.time()
data = CreateSeuratObject(counts = counts, project = project_name, min.cells = 3, min.features = 200)
name = paste("Seurat", project_name, 'CreateSeuratObject', sep=".")
track_mem(mem_info_start, file.path(write_dir, name))
time_dict['CreateSeuratObject'] = Sys.time() - t_ref
mem_info_start = gc(reset=TRUE)

t_ref = Sys.time()
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
name = paste("Seurat", project_name, 'PercentageFeatureSet', sep=".")
track_mem(mem_info_start, file.path(write_dir, name))
time_dict['PercentageFeatureSet'] = Sys.time() - t_ref
mem_info_start = gc(reset=TRUE)


t_ref = Sys.time()
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
name = paste("Seurat", project_name, 'NormalizeData', sep=".")
track_mem(mem_info_start, file.path(write_dir, name))
time_dict['NormalizeData'] = Sys.time() - t_ref
mem_info_start = gc(reset=TRUE)


t_ref = Sys.time()
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(data), 10)
name = paste("Seurat", project_name, 'FindVariableFeatures', sep=".")
track_mem(mem_info_start, file.path(write_dir, name))
time_dict['FindVariableFeatures'] = Sys.time() - t_ref
mem_info_start = gc(reset=TRUE)


t_ref = Sys.time()
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
name = paste("Seurat", project_name, 'ScaleData', sep=".")
track_mem(mem_info_start, file.path(write_dir, name))
time_dict['ScaleData'] = Sys.time() - t_ref
mem_info_start = gc(reset=TRUE)


t_ref = Sys.time()
data <- RunPCA(data, features = VariableFeatures(object = data))
name = paste("Seurat", project_name, 'RunPCA', sep=".")
track_mem(mem_info_start, file.path(write_dir, name))
time_dict['RunPCA'] = Sys.time() - t_ref
mem_info_start = gc(reset=TRUE)


t_ref = Sys.time()
data <- FindNeighbors(data, dims = 1:10)
name = paste("Seurat", project_name, 'FindNeighbors', sep=".")
track_mem(mem_info_start, file.path(write_dir, name))
time_dict['FindNeighbors'] = Sys.time() - t_ref
mem_info_start = gc(reset=TRUE)


t_ref = Sys.time()
data <- FindClusters(data, resolution = 0.5)
name = paste("Seurat", project_name, 'FindClusters', sep=".")
track_mem(mem_info_start, file.path(write_dir, name))
time_dict['FindClusters'] = Sys.time() - t_ref
mem_info_start = gc(reset=TRUE)


t_ref = Sys.time()
data <- RunUMAP(data, dims = 1:10)
name = paste("Seurat", project_name, 'RunUMAP', sep=".")
track_mem(mem_info_start, file.path(write_dir, name))
time_dict['RunUMAP'] = Sys.time() - t_ref
mem_info_start = gc(reset=TRUE)


# -- find cluster-specific markers: ------------------------------------------------------
# ---- find markers for every cluster compared to all remaining cells, 
# ---- report only positive markers
t_ref = Sys.time()
data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
name = paste("Seurat", project_name, 'FindAllMarkers', sep=".")
track_mem(mem_info_start, file.path(write_dir, name))
time_dict['FindAllMarkers'] = Sys.time() - t_ref
mem_info_start = gc(reset=TRUE)


# -----------------------------------------------------------------------------------------------
name = paste("Seurat", project_name, 'TotalMemory', sep=".")
track_mem(mem_info_start, file.path(write_dir, name))
time_dict['tf'] = Sys.time()
time_dict['total_diff'] = time_dict['tf'] - time_dict['t0']


fp = file.path(write_dir, paste("Seurat", project_name, 'TotalTime', 'csv', sep="."))
write.csv(as.data.frame(time_dict), fp)

