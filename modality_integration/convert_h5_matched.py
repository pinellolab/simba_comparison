import scanpy as sc
import pandas as pd
import sys

input_h5ad = sys.argv[1]
output_file = sys.argv[2]

adata = sc.read_h5ad(input_h5ad)
df = pd.DataFrame(adata.X, index = adata.obs_names)
df.rename(lambda x: "PC" + str(x) , axis = 'columns', inplace = True)
df['cell_type'] = adata.obs.celltype
df['batch'] = adata.obs.entity_group
df.to_csv(output_file, sep = "\t")
