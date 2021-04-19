import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import time
from datetime import timedelta
import scanpy as sc
import os
import sys
import random

from sklearn.metrics import silhouette_score

# asw_batch_norm_sub = 1 - asw_batch_norm
def silhouette_coeff_ASW(adata, method_use='raw',save_dir='', save_fn='', percent_extract=0.8, metric = "euclidean"):
    random.seed(0)
    asw_b = []
    asw_b_sub = []
    asw_ct = [] 
    iters = []
    for i in range(20):
        iters.append('iteration_'+str(i+1))
        rand_cidx = np.random.choice(adata.obs_names, size=int(len(adata.obs_names) * percent_extract), replace=False)
#         print('nb extracted cells: ',len(rand_cidx))
        adata_ext = adata[rand_cidx,:]
        asw_batch = silhouette_score(adata_ext.X, adata_ext.obs['batch'], metric = metric)
        asw_celltype = silhouette_score(adata_ext.X, adata_ext.obs['cell_type'], metric = metric)
        
        asw_b.append(asw_batch)
        asw_b_sub.append(1-asw_batch)
        asw_ct.append(asw_celltype)
    
#     iters.append('median_value')
#     asw_fscore.append(np.round(np.median(fscoreASW),3))
#     asw_bn.append(np.round(np.median(asw_batch_norm),3))
#     asw_bn_sub.append(np.round(1 - np.median(asw_batch_norm),3))
#     asw_ctn.append(np.round(np.median(asw_celltype_norm),3))
    df = pd.DataFrame({'asw_batch':asw_b, 'asw_batch_sub': asw_b_sub,
                       'asw_celltype': asw_ct,
                       'method_use':np.repeat(method_use, len(asw_ct))})
    df = df.append(df.median(axis = 0), ignore_index=True)
    df.loc[max(df.index), "method_use"] = "median"

    df.to_csv(save_dir + save_fn + '.txt', sep = "\t")
    print('Save output of pca in: ',save_dir)
    print(df.values.shape)
    print(df.keys())
    return df
        
def createAnnData(myDatafn):
    
    myData = pd.read_table(myDatafn, index_col = 0)
    bex = ['batch','Batch','Batchlb','batchlb','BATCH']
    ib = np.isin(myData.keys(), bex)
    cex = ['celltype','CellType','cell_type','Cell_Type','ct']
    ict = np.isin(myData.keys(), cex)
    adata = sc.AnnData(myData.values[:,0:np.shape(myData.values)[1]-2])
    adata.obs_names = myData.index
    adata.obs['batch'] = myData.values[:, np.where(ib)[0][0]]  # factor function in R
    adata.obs['cell_type'] = myData.values[:, np.where(ict)[0][0]]
    print(adata)
    return adata

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
#print(sc.logging.print_versions())

method_use = sys.argv[1]
pca_file = sys.argv[2]
save_dir = sys.argv[3]
emb_type = sys.argv[4]
dissim_metric = sys.argv[5]
print("Saving to ", save_dir)

#if not os.path.exists(save_dir+'/ASW/'): os.makedirs(os.path.join(save_dir,'/ASW/')) 

# _pca.csv in data_dir
#fls = [ f for f in os.listdir(data_dir) if f.endswith("_pca.csv") & os.path.isfile(os.path.join(data_dir,f)) ]
#fls

def main(f = pca_file):
    #final_ls = []
    print('Extract asw for ', method_use)
    save_fn = method_use + '_ASW_' + emb_type + '_' + dissim_metric
    adata = createAnnData(f)
    asw_val = silhouette_coeff_ASW(adata, method_use, save_dir, 
                                          save_fn, percent_extract=0.8, 
                                          metric = dissim_metric)
    #final_ls.append(asw_val)
    

#result = pd.concat(final_ls)   
#result.to_tsv(os.path.join(save_dir, method_use + 'ASW.txt'))

main(pca_file)
