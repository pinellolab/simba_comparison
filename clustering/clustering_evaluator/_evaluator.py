
from sklearn.metrics.cluster import (
    adjusted_rand_score,
    adjusted_mutual_info_score,
    homogeneity_score,
)

import pandas as pd


class Evaluator:
    def __init__(
        self, adata, ref_key="celltype", cluster_keys=["louvain", "kmeans", "h_clust"]
    ):

        self.adata = adata
        self._ref_key = ref_key
        self._cluster_keys = cluster_keys

        self.scores = {"ARI": {}, "AMI": {}, "Homogeneity": {}}

    def ARI(self):
        """adjusted rank index score"""
        for clust_key in self._cluster_keys:
            self.scores["ARI"][clust_key] = adjusted_rand_score(
                self.adata.obs[self._ref_key], self.adata.obs[clust_key]
            )

    def AMI(self):
        """adjusted mutual information score"""
        for clust_key in self._cluster_keys:
            self.scores["AMI"][clust_key] = adjusted_mutual_info_score(
                self.adata.obs[self._ref_key],
                self.adata.obs[clust_key],
                average_method="arithmetic",
            )

    def Homogeneity(self):
        """homogeneity score"""
        for clust_key in self._cluster_keys:
            self.scores["Homogeneity"][clust_key] = homogeneity_score(
                self.adata.obs[self._ref_key], self.adata.obs[clust_key]
            )

    def __call__(self):
        self.ARI()
        self.AMI()
        self.Homogeneity()
        self.scores = pd.DataFrame(self.scores)
        self.adata.uns["clustering_scores"] = self.scores