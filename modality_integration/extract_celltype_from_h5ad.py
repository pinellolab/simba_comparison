import scanpy as sc
import pandas as pd
import sys

input_h5ad = "/data/pinello/PROJECTS/2019_08_Embedding/SIMBA_RESULTS/multiome_10xpbmc_10k/result_10xpbmc10k_integration_v1/adata_all.h5ad"
output_file = "data/pbmc/cell_type_h5ad.txt"

print("convert h5ad ", input_h5ad, " into ", output_file)

adata = sc.read_h5ad(input_h5ad)

df = adata.obs.celltype.to_frame("cluster")
df.reset_index(inplace = True)
df = df.rename(columns = {'index' : 'bc'})
df['bc'] = list(map(lambda x: str(x).split("_")[0], df['bc']))


df.to_csv(output_file, sep = "\t", header = True, index = False)
