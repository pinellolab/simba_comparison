import bbknn
import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
import os
import sys
import time
from datetime import timedelta

sc.settings.verbosity = 3

dirname = os.getcwd()
print(dirname)
data_dir = os.path.join(dirname, './data/')
# if not os.path.exists('./results_cellatlas/'): os.makedirs('./results_cellatlas/')

if not os.path.exists('./results_cellatlas/bbknn_results_v2/'): os.makedirs('./results_cellatlas/bbknn_results_v2/')

expr_data_path = sys.argv[1]
metadata_path = sys.argv[2]
output_prefix = sys.argv[3]
save_dir = sys.argv[4]

def save_images(basename):
    if not os.path.exists(save_dir): os.makedirs(save_dir)
    pl.savefig(save_dir + "/" + outname + ".png" , dpi=150)
    pl.close()

adata1 = sc.read_h5ad(os.path.join(data_dir,'hvg_dataset2_cellatlas.h5ad')) # filtered data, keep only hvg genes 6954 × 1328 

print(adata1)
# sc.pp.log1p(adata1)
# sc.pp.scale(adata1)
sc.tl.leiden(adata1)
bbknn.ridge_regression(adata1, batch_key = "batch", confounder_key = "leiden")
sc.tl.pca(adata1, svd_solver='arpack')
adata1.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
sc.pl.pca_variance_ratio(adata1, log=True)
# sc.pp.neighbors(adata1,n_neighbors=10, n_pcs=20)

sc.pp.neighbors(adata1,n_pcs=num_pcs, n_neighbors=20)
sc.tl.umap(adata1)
sc.tl.louvain(adata1)
color_group = ['cell_type','batchlb','louvain']
sc.pl.umap(adata1, color=color_group, save = )

del adata1.uns['louvain_colors']

## Options
trim = 50 # Trim the neighbors of each cell to these many top connectivities
neighbors_within_batch = 5 # How many top neighbors to report for each batch; total # of neighbors = this * #batches
num_pcs = 20
approx = True # annoy's approximate neighbor finding. Quicker run time for large dataset
use_faiss = True # when approx = False and metric = "euclidean", use faiss package to compute nearest neighbors. Improve numerical precision

t3 = time.time()
adata_bbknn_trim = bbknn.bbknn(adata1,copy=True,neighbors_within_batch=neighbors_within_batch,trim=trim,n_pcs=num_pcs,batch_key='batch') #approx=False,
t4 = time.time()
print('Took '+str(timedelta(seconds=t4-t3)))

sc.tl.pca(adata_bbknn_trim, svd_solver='arpack',n_comps=20)
adata_bbknn_trim.obsm['X_pca'] *= -1

def write_to_csv(mat, genesname, cellsname, filename, save_dir):
    if isinstance(mat, np.ndarray):
        df = pd.DataFrame(mat, columns=genesname, index=cellsname)
    else:
        df = pd.DataFrame(mat.toarray(), columns=genesname, index=cellsname)        
    
    df.to_csv(save_dir+filename)  
    
filename = 'murine-atlas_BBKNN_pca.txt'
coln_pca = []
for i in range(adata_bbknn.obsm['X_pca'][:,:20].shape[1]):
    coln_pca.append("X_pca"+str(i+1))
    

write_to_csv(adata_bbknn_trim.obsm['X_pca'], coln_pca, adata_bbknn_trim.obs_names,filename, save_dir)

