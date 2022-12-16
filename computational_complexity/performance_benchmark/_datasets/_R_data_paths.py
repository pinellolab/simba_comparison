
# -- import packages: --------------------------------------------------------------------
import os


# -- supporting functions: ---------------------------------------------------------------
def _format_data_paths(self, path):
    """Format paths from DataConfig"""
    for key in ["h5ad", "obs", "var"]:
        if key in ["obs", "var"]:
            path_ = path + ".{}_df.csv".format(key)
        else:
            path_ = path + ".h5ad"
        setattr(self, "{}_path".format(key), path_)


# -- API-facing class: -------------------------------------------------------------------
class RDataPaths:
    """Class for managing paths to R data derived from .h5ad"""
    def __init__(self,
                 name="unspecified",
                 path="/home/mvinyard/SIMBA/data/",
                 n_downsample=0,
                ):
        
        _format_data_paths(self, path)
        self.N = n_downsample
        self.name = name
        
    def verify_paths(self):
        """Verify that created paths exist for adata components."""
        self._verified = all(
            [
                os.path.exists(self.h5ad_path),
                os.path.exists(self.obs_path),
                os.path.exists(self.var_path),
            ]
        )
        return self._verified