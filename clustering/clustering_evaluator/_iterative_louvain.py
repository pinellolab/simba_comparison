
import numpy as np
import scanpy as sc
import pandas as pd

class DynamicClusteringResolution:
    def __parse__(self, kwargs, ignore=["self"]):
        for k, v in kwargs.items():
            if not k in ignore:
                if k == "scan":
                    v = list(v)
                setattr(self, "_{}".format(k), v)

    def __init__(self, target, scan=(1, 5)):
        self.__parse__(locals())
        self.current = np.mean(self._scan)

    def _update(self, predicted):
        if predicted > self._target:
            self._scan[1] = self.current
        else:
            self._scan[0] = self.current
        self.current = np.mean(self._scan)

    def evaluate(self, predicted):
        if predicted != self._target:
            self._update(predicted)


class IterativeLouvain:
    def __init__(self, adata, target=9, scan=(0.1, 2)):
        self.adata = adata
        self._scan = scan
        self._target = target
        self.DCR = DynamicClusteringResolution(target=target, scan=scan)
        self.RunInfo = {}

    def _has_neighbors(self, adata):
        return "neighbors" in adata.uns_keys()

    def _neighbors(self, adata, n_neighbors=15, use_rep="X"):
        if not self._has_neighbors(adata):
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep)

    def _log(self):
        pass

    def _forward(self, resolution):
        self._neighbors(self.adata, n_neighbors=15, use_rep="X")
        sc.tl.louvain(self.adata, resolution=resolution)
        n_clusters = self.adata.obs["louvain"].nunique()
        self.DCR.evaluate(predicted=n_clusters)
        self.run_info = {"resolution": resolution, "n_clusters": n_clusters}

        return n_clusters

    def run(self, max_iters=20):

        for i in range(max_iters):
            n_pred = self._forward(resolution=self.DCR.current)
            self.RunInfo[i] = self.run_info
            if n_pred == self._target:
                break

        self.RunInfo = pd.DataFrame.from_dict(self.RunInfo).T
                
