
import simba as si

def rna_10x_pbmcs_3k(quiet=False):
    adata = si.datasets.rna_10xpmbc3k()
    if not quiet:
        print(adata)
    return adata
