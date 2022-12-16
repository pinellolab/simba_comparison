

# -- import local dependencies: ----------------------------------------------------------
from .._tools._seurat_resource_collector import SeuratResourceCollector
from .._tools._read_pickle_file import read_pkl


# -- import packages: --------------------------------------------------------------------
from abc import ABC, abstractmethod
import pandas as pd
import numpy as np
import glob
import os


# -- supporting functions: ---------------------------------------------------------------
def is_uniform(x):
    return len(set(x)) == 1


# -- base class: -------------------------------------------------------------------------
class Results(ABC):
#     def _augment_fnames(self, fnames):
#         fname_lens = [len(fname) for fname in fnames]
#         if not is_uniform(fname_lens):
#             min_idx = np.argmin(fname_lens)
#             fnames[min_idx] = fnames[min_idx] + "_0"
#         return fnames

    @abstractmethod
    def mk_filename(self, path) -> str:
        return fname

    def compose_file_df(self):
        
        self.results_paths = glob.glob(self.path)
        fnames = [self.mk_filename(path) for path in self.results_paths]

#         fnames = self._augment_fnames(fnames)
        f_iters = [int(fname.split(".")[self.iter_pos]) for fname in fnames]
        self.file_df = (
            pd.DataFrame(
                {"iter": f_iters, "fnames": fnames, "paths": self.results_paths}
            )
            .set_index("iter")
            .sort_index()
        )

    @abstractmethod
    def read_pkl(self):
        pass

    def compose_FileDict(self):
        self.FileDict = {
            i: self.read_pkl(self.file_df["paths"][i]) for i in self.file_df.index
        }

    def __call__(self):
        self.compose_file_df()
        self.compose_FileDict()
        return self.FileDict


# -- detailed sub-classes: ---------------------------------------------------------------
class SIMBAResults(Results):
    def __init__(self, path):
        self.path = path
        self.iter_pos = -2

    def mk_filename(self, path):
        return os.path.basename(path)

    def read_pkl(self, path):
        return read_pkl(path)


class ScanpyResults(Results):
    def __init__(self, path):
        self.path = path
        self.iter_pos = -2
        
    def mk_filename(self, path):
        return os.path.basename(path)

    def read_pkl(self, path):
        return read_pkl(path)

class SeuratResults(Results):
    def __init__(self, path):
        self.path = path
        self.iter_pos = -1

    def mk_filename(self, path):
        return os.path.basename(os.path.dirname(path))

    def read_pkl(self, path):
        return SeuratResourceCollector(path)