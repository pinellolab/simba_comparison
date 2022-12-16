

# -- import packages: -------------------------------------------------------
import os
from tqdm.notebook import tqdm


# -- import local dependencies: ---------------------------------------------
from . import _scanpy_funcs as funcs
from ._format_scanpy_results import format_scanpy_results
from ... import _tools as tl



def _n_cells(adata):
    return adata.shape[0]


def _run_scanpy_traced(adata):
    
    adata_ = adata.copy()

    run_scanpy = funcs.TracedScanpy()
    run_scanpy(adata=adata_)
    resource_df = run_scanpy.outs[1]
    
    return run_scanpy, resource_df


# -- API-facing exe and collection of results: ------------------------------
def run_scanpy(adata,
               n_iters=5,
               method_dir="./scanpy_benchmark_results",
               dataset="unspecified_dataset",
               n_cells=None,
               parent_path=os.getcwd(),
               debug=False,
             ):
    
    outs = {}
    
    if not n_cells:
        n_cells = _n_cells(adata)
        
    for i in tqdm(range(n_iters)):
        filename = tl.mk_save_path(method_dir, dataset, ext=".zip")
        run_outs = _run_scanpy_traced(adata)
        if debug:
            return run_outs
        outs[i] = format_scanpy_results(run = run_outs,
                              save_filename = filename,
                              dataset = dataset,
                              n_cells = n_cells,
                              run_iter = i,
                              return_saved_outs=True
                             )
    
    return outs
    