import anndata as a
import glob
import os


def fetch_name(path):
    return ".".join(os.path.basename(path).split(".")[2:4])


def report_datasets(DataDict):
    for name, adata in DataDict.items():
        print(
            "{:<30} | celltypes: {:<2} | dataset: {:<14} | features: {}".format(
                adata.uns["name"],
                adata.uns["n_clusters"],
                adata.uns["dataset"],
                adata.uns["features"],
            )
        )


def load_data(path, silent=False):

    """
    Parameters:
    -----------
    path

    Notes:
    ------
    Celltypes are encoded as adata.obs['celltype']
    """

    h5ad_paths = glob.glob(path + "*.h5ad")

    DataDict = {}

    for path in h5ad_paths:
        name = fetch_name(path)
        adata = a.read_h5ad(path)
        adata.uns["n_clusters"] = adata.obs["celltype"].nunique()
        adata.uns["name"] = name
        adata.uns["path"] = path
        adata.uns["dataset"] = name.split(".")[0]
        adata.uns["features"] = name.split(".")[1:][0].split("_")
        DataDict[name] = adata

    if not silent:
        report_datasets(DataDict)

    return DataDict