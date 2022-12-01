

from ._clustering import Clustering
from ._evaluator import Evaluator
import scanpy as sc
from tqdm import tqdm

class ClusteringEvaluator:
    def __init__(self, adata, n_clusters=9, random_state=2022):
        self.adata = adata
        self.n_clusters = n_clusters
        self.random_state = random_state

    def run(self):

        self.clustering = Clustering(
            self.adata, n_clusters=self.n_clusters, random_state=self.random_state
        )
        self.clustering()
        self.evaluator = Evaluator(self.adata)
        self.evaluator()

        self.scores = self.evaluator.scores

    def quick_plot(self):
        sc.pl.umap(self.adata, color=["louvain", "kmeans", "h_clust"], frameon=False)
        

def evaluate_embedding(DataDict, n_clusters=None, random_state=2022, plot=False):
    
    for name, adata in tqdm(DataDict.items(), ncols=50):
        ce = ClusteringEvaluator(adata, n_clusters, random_state=random_state)
        ce.run()
    
        if plot:
            ce.quick_plot()