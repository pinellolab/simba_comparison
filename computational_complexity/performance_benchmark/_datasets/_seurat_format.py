

from anndata import AnnData


class SeuratFormat:
    def __init__(self, adata, save=False):

        self.adata = adata
        self._save = save
        self.var_df = self.adata.var.copy()
        self.obs_df = self.adata.obs.copy()

    def save_paths(self):
        if self._save:
            self.SavePaths = {
                "h5ad": ".".join([self._save, "h5ad"]),
                "var": ".".join([self._save, "var_df.csv"]),
                "obs": ".".join([self._save, "obs_df.csv"]),
            }

    def format_adata(self):
        adata_ = AnnData(self.adata.X)

        if self._save:
            adata_.write_h5ad(self.SavePaths["h5ad"])
            print("\nAnnData written to: {}".format(self.SavePaths["h5ad"]))

    def format_var_df(self):

        self.var_df = self.var_df.reset_index().rename({"index": "gene_names"}, axis=1)

        if self._save:
            self.var_df.to_csv(self.SavePaths["var"])
            print(".var DataFrame written to: {}".format(self.SavePaths["var"]))

    def format_obs_df(self):

        self.obs_df = self.obs_df.reset_index().rename({"index": "barcodes"}, axis=1)

        if self._save:
            self.obs_df.to_csv(self.SavePaths["obs"])
            print(".obs DataFrame written to: {}".format(self.SavePaths["obs"]))

    def __call__(self):

        self.save_paths()
        self.format_adata()
        self.format_var_df()
        self.format_obs_df()