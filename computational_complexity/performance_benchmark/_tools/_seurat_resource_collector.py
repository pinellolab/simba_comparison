


# -- import packages: --------------------------------------------------------------------
import os
import glob
import pandas as pd


# -- supporting functions: --------------------------------------------------------
def bits_to_bytes(n_bits):
    return (n_bits * 8) / (1024 * 1024)


# -- main module class: ------------------------------------------------------------------
class SeuratResourceCollector:
    def __glob__(self, write_dir):
        return glob.glob(write_dir + "/*")

    def __parse__(self):
        out_dict = {}
        for path in self.outs:
            name = os.path.basename(path).split(".")[-2]
            if not "Time" in name:
                out_dict[name] = pd.read_csv(path, index_col=0)
            else:
                time = pd.read_csv(path, index_col=0)
        mem_df = pd.concat(list(out_dict.values()), axis=1)
        mem_df.columns = list(out_dict.keys())

        self.mem = mem_df.T
        self.time = time

    def __init__(self, write_dir):
        self.outs = self.__glob__(write_dir)
        self.__parse__()
        
    def compute_mb_used(self):

        return (
            (self.mem["mem_end_max"] - self.mem["mem_start_max"])
            .apply(bits_to_bytes)
#             .iloc[1:]
        )