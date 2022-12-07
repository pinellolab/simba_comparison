

# -- import packages: --------------------------------------------------------------------
import pandas
import anndata
import pandas as pd
import numpy as np


# -- supporting function: ----------------------------------------------------------------
def _format_clustering_scores(adata: anndata.AnnData, name: str):

    """
    adata

    name: dictionary key
    """

    tmp_df1 = pd.DataFrame(name.split("."), index=["dataset", "features"]).T

    tmp_df2 = adata.uns["clustering_scores"].unstack().to_frame().reset_index()
    idx = tmp_df2["level_0"].astype(str) + "_" + tmp_df2["level_1"].astype(str)
    score_df = pd.DataFrame(tmp_df2[0].values, index=idx).T
    return pd.concat([tmp_df1, score_df], axis=1)


# -- main module function: ---------------------------------------------------------------
def format_clustering_scores(DataDict: dict) -> pandas.DataFrame:

    """
    Parameters:
    -----------
    DataDict
        Dictionary of adata objects
        type: dict

    Returns:
    --------
    score_df
        type: pandas.DataFrame

    Notes:
    ------

    """

    data_keys = np.sort(list(DataDict.keys())).tolist()
    scores = [_format_clustering_scores(DataDict[name], name) for i, name in enumerate(data_keys)]
    return pd.concat(scores).reset_index(drop=True)