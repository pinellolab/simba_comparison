
from sklearn.cluster import KMeans, AgglomerativeClustering
from ._iterative_louvain import IterativeLouvain
import pandas as pd


class Clustering:
    def __init__(self, adata, n_clusters=None, random_state=2022):
        self.adata = adata
        
        if not n_clusters:
            n_clusters = adata.uns['n_clusters']
        
        self.n_clusters = n_clusters
        self.random_state = random_state

    def Louvain(self):
        self._louvain = IterativeLouvain(self.adata, target=self.n_clusters, scan=(0.1, 3))
        self._louvain.run()

    def KMeans(self):
        self._k_means = KMeans(n_clusters=self.n_clusters, random_state=self.random_state).fit(
            self.adata.X
        )
        labels = self._k_means.labels_
        idx = self.adata.obs.index
        self.adata.obs["kmeans"] = pd.Series(labels, index=idx).astype("category")

    def Hierarchical(self):
        self._h_clust = AgglomerativeClustering(n_clusters=self.n_clusters).fit(
            self.adata.X
        )
        labels = self._h_clust.labels_
        idx = self.adata.obs.index
        self.adata.obs["h_clust"] = pd.Series(labels, index=idx).astype("category")

    def __call__(self):
        
        self.Louvain()
        self.KMeans()
        self.Hierarchical()