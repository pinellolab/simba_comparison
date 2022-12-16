
import os

from ._seurat_format import SeuratFormat
from ._R_data_paths import RDataPaths


def downsample(adata, n, random_state=None, quiet=False, base_path_for_seurat=False, name=False, **kwargs):

    """
    Use the pandas `sample` function to downsample cells from adata.

    Parameters:
    -----------
    adata

    n

    random_state
    
    base_path_for_seurat
        if not False (default), str should serve as the base path.
        E.g.: base_path_for_seurat= './result_simba/data/rna_mouse_neurons'

    **kwargs

    Returns:
    --------
    adata
        Downsampled AnnData object of shape: n x n_genes
    """

    if n == adata.obs.shape[0]:
        return adata
    
    
    idx = adata.obs.sample(n=n, random_state=random_state, **kwargs).index
    
    adata_ = adata[idx].copy()
    
    if not quiet:
        print(adata_)
        
    if base_path_for_seurat:
        path = "{}.{}cells".format(base_path_for_seurat, n)
        seurat_formatter = SeuratFormat(
            adata_, save=path
        )
        seurat_formatter()
        
        if not name:
            name = os.path.basename(base_path_for_seurat)
        return RDataPaths(name=name, path=path)
    else:    
        return adata_