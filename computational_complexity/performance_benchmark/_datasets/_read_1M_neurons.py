
import os
import simba as si
import pandas as pd

def read_1M_neurons(data_dir = "../data/1M_neurons", quiet=False):
    
    h5_path = "1M_neurons_filtered_gene_bc_matrices_h5.h5"
    clusters_path = "clustering/graphclust/clusters.csv"
    
    adata = si.read_10x_h5(os.path.join(data_dir, h5_path))
    adata.var_names_make_unique()
    
    df_metadata = pd.read_csv(os.path.join(data_dir, clusters_path), index_col=0)
    df_metadata.index.name = None
    df_metadata['Cluster'] = 'cluster_' + df_metadata['Cluster'].astype(str)
    adata.obs = df_metadata
    
    if not quiet:
        print(adata)
    
    return adata