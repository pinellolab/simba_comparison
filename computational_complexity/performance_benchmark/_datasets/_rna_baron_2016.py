
import simba as si

def rna_baron_2016(quiet=False):
    adata = si.datasets.rna_baron2016()
    if not quiet:
        print(adata)
    return adata