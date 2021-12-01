import scanpy as sc
import pandas as pd
import sys

input_h5ad = sys.argv[1]
output_file = sys.argv[2]
umap_file = sys.argv[3]

print("convert h5ad ", input_h5ad, " into ", output_file)

adata = sc.read_h5ad(input_h5ad)
df = pd.DataFrame(adata.X, index = adata.obs_names)
df.rename(lambda x: "D" + str(x) , axis = 'columns', inplace = True)
df['cell_type'] = adata.obs.celltype
df['cell_type'] = list(map(lambda x: str(x).replace(" ", "_"), df['cell_type']))
df['batch'] = adata.obs.entity_group
df.to_csv(output_file, sep = "\t")

df_umap = pd.DataFrame(adata.obsm["X_umap"], index = adata.obs_names)
df_umap.columns = ["UMAP_1", "UMAP_2"]
df_umap['cell_type'] = adata.obs.celltype
df_umap['cell_type'] = list(map(lambda x: str(x).replace(" ", "_"), df_umap['cell_type']))
df_umap['batch'] = adata.obs.entity_group
df_umap.to_csv(umap_file, sep = "\t")

