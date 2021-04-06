import scanpy as sc
import pandas as pd
import sys

input_h5ad = sys.argv[1]
output_path = sys.argv[4]

adata = sc.read_h5ad(input_h5ad)

df = pd.DataFrame(adata.obsm.get("X_umap"), columns = ["UMAP_1", "UMAP_2"], index = adata.obs_names)
df['cell_type'] = adata.obs.celltype
df['cell_type'] = list(map(lambda x: str(x).replace(" ", "_"), df['cell_type']))
df['modality'] = adata.obs.entity_group
df.to_csv(output_path, sep = "\t")
