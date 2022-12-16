
#
#
#
#

# -- import packages: --------------------------------------------------------------------
import os
import subprocess


# -- import local dependencies: ----------------------------------------------------------
from ..._tools._mk_save_path import mk_save_path
from ..._tools._seurat_resource_collector import SeuratResourceCollector

# -- Seurat controller class: ------------------------------------------------------------
class Seurat:
    def __parse__(self, kwargs, ignore=["self"], hide=["dataset", "config"]):

        self._kwargs = {}
        for k, v in kwargs.items():
            if not k in ignore:
                if k in hide:
                    self._kwargs["_{}".format(k)] = v
                    setattr(self, "_{}".format(k), v)
                else:
                    self._kwargs[k] = v
                    setattr(self, k, v)

    def __init__(self,
                 config={},
                 R="/home/mvinyard/.anaconda3/envs/seurat/bin/Rscript",
                ):

        self._R = R
        self._N = config.N
        self.dataset = config.name
        _seurat_path = os.path.join(os.path.dirname(__file__), "_run_seurat.R")
        
        for obj in ["h5ad", "obs", "var"]:
            setattr(self, "_{}_path".format(obj), getattr(config, "{}_path".format(obj)))
                
        self.__parse__(locals())

        self._write_dir = mk_save_path(
            method_dir="seurat_benchmark_results",
            dataset_dir="seurat.{}".format(self.dataset),
            ext="",
        )

    def configure(self):
        self.command = command = [
            self._R,
            self._seurat_path,
            self._h5ad_path,
            self._obs_path,
            self._var_path,
            self._write_dir,
            self.dataset,
            str(self._N),
        ]

    def forward(self):
        self.return_code = subprocess.run(
            self.command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        self.stdout = self.return_code.stdout.decode("utf-8").split("\n")

        self.stderr = self.return_code.stderr.decode("utf-8").split("\n")

    def run(self):
        self.configure()
        self.forward()
        
        
def _run_seurat(config,
                n_iters: int = 4,
                R: str = "/home/mvinyard/.anaconda3/envs/seurat/bin/Rscript",
                notebook: bool = True,
               ):    
    """
    API-facing function to run Seurat
    
    config_dict
        This is wher you specify paths to input data.
    """
    
    if notebook:
        from tqdm.notebook import tqdm
    else:
        from tqdm import tqdm

    seurat_outs = {}
    for i in tqdm(range(n_iters)):
        seurat = Seurat(config, R=R)
        seurat.run()
        try:
            seurat_outs[i] = SeuratResourceCollector(seurat._write_dir)
        except:
            print(" - [ WARNING ] | There was a problem writing / collecting Seurat iteration: {}...".format(i))
        
    return seurat_outs, seurat
