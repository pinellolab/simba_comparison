
# -- import packages: --------------------------------------------------------------------
import os
import pickle
import pandas as pd
import licorice_font
from anndata import AnnData


# -- import local dependencies: ----------------------------------------------------------


# -- supporting functions: ---------------------------------------------------------------
import pandas as pd
import glob
import os

import datetime

def _from_iso(x):
    s = datetime.datetime.fromisoformat(x).second / (60*60)
    m = datetime.datetime.fromisoformat(x).minute / 60
    h = datetime.datetime.fromisoformat(x).hour
    return h + m + s

def load_memory_trace(path=None):
    
    """
    By default, will read the most recent. 
    """
    
    if not path:
        path="simba.background_memory_trace*.txt"
    
        memory_trace_files = {}

        for p in glob.glob(path):
            split_p = os.path.basename(p).split(".")
            if not len(split_p) == 4:
                version = 0
            else:
                version = split_p[-2]
            memory_trace_files[int(version)] = p

        mem_paths = pd.DataFrame.from_dict(memory_trace_files, orient="index").reset_index()
        mem_paths.columns = ["iter", "path"]
        mem_paths['iter'] = mem_paths['iter'].astype(int)

        mem_path_df = mem_paths.sort_values('iter').reset_index(drop=True)

        use_path = mem_path_df.iloc[mem_path_df['iter'].iloc[-1]]['path']
        
    else:
        use_path = path
        
    mem_df = pd.read_csv(use_path, usecols=range(12))
    
    return mem_df, use_path

def rm_adata_simba(runner):
    for attr in ["_adata", "_adata_CG"]:
        if hasattr(runner, attr):
            delattr(runner, attr)


def rm_adata_scanpy(runner):
    for attr in ["_adata"]:
        if hasattr(runner, attr):
            delattr(runner, attr)

    runner.outs = list(runner.outs)
    for obj in runner.outs:
        if isinstance(obj, AnnData):
            del obj
    return runner


def to_pickle(obj, filename, protocol=pickle.HIGHEST_PROTOCOL, verbose=True):
    pickle.dump(obj, open(filename, "wb"), protocol=protocol)
    if verbose:
        msg = licorice_font.font_format("Successfully saved to", ["BOLD"])
        print("{}: {}".format(msg, filename))


def size_of(obj: dict):
    """Returns object size in bytes"""
    s = 0
    size_dict = {}
    for k, v in obj.items():
        size_dict[k] = _s = os.sys.getsizeof(v)
        s += _s
    return s, size_dict

def format_graph_time(gen_graph):

    graph_total_t = gen_graph.resource_df["t"].values[0]
    edge_write_t = gen_graph.outs["edge_write_time"]
    graph_gen_t = graph_total_t - edge_write_t

    return pd.DataFrame(
        [edge_write_t, graph_gen_t], index=["write_graph", "gen_graph"], columns=["t"]
    )

# -- manager class: ----------------------------------------------------------------------
class FormattedResults:
    def __init__(self, outs):
        self.outs = outs

    # -- SIMBA: --------------------------------------------------------------------------
    # ((runner, resource_df), mem), t = outs
    @property
    def simba_runner(self):
        _runner = self.outs[0][0][0]
        rm_adata_simba(_runner)
        return _runner
    
    @property
    def simba_landmark_datetimes(self):
        return self.outs[0][0][2]
        
    @property
    def simba_memory_trace_df(self):
        return load_memory_trace()

    @property
    def simba_resource_df(self):
        _resource_df = self.outs[0][0][1]
        _graph_gen_df = format_graph_time(self.simba_runner._gen_graph)
        return pd.concat([_resource_df, _graph_gen_df], axis=0)

    @property
    def simba_total_mem(self):
        return self.outs[0][1].components

    @property
    def simba_total_time(self):
        return self.outs[1]

    # -- Scanpy: -------------------------------------------------------------------------
    @property
    def scanpy_resource_df(self):
        return self.outs[1]

    @property
    def scanpy_runner(self):
        return rm_adata_scanpy(self.outs[0])

    @property
    def scanpy_total_mem(self):
        return self.scanpy_runner.m.components

    @property
    def scanpy_total_time(self):
        return self.scanpy_runner.t

    def metadata(self, dataset=None, n_cells=None, run_iter=None):
        return {"dataset": dataset, "n_cells": n_cells, "run_iter": run_iter}

    # -- do formatting: ------------------------------------------------------------------
    def format_simba_results(
        self, dataset=None, n_cells=None, run_iter=None, include_runner=False
    ):
        
        mem_trace_df, use_path = self.simba_memory_trace_df
        landmarks = self.simba_landmark_datetimes
        
        mem_trace_df = mem_trace_df.loc[
            (mem_trace_df["datetime"] > landmarks["t0"])
            & (mem_trace_df["datetime"] < landmarks["tf"])
        ]

        self.formatted_results = {
            "mem_total": self.simba_total_mem,
            "time_total": self.simba_total_time,
            "resource_df": self.simba_resource_df,
            "landmark_datetimes": landmarks,
            "memory_trace_df": mem_trace_df,
            "metadata": self.metadata(
                dataset=dataset,
                n_cells=n_cells,
                run_iter=run_iter,
            ),
        }

        if include_runner:
            self.formatted_results["runner"] = self.simba_runner
            
    def format_scanpy_results(
        self, dataset=None, n_cells=None, run_iter=None, include_runner=False
    ):

        self.formatted_results = {
            "mem_total": self.scanpy_total_mem,
            "time_total": self.scanpy_total_time,
            "resource_df": self.scanpy_resource_df,
            "metadata": self.metadata(
                dataset=dataset,
                n_cells=n_cells,
                run_iter=run_iter,
            ),
        }

        if include_runner:
            self.formatted_results["runner"] = self.scanpy_runner

    # -- save formatted results ----------------------------------------------------------
    def save(self, filename):
        self.filename = filename
        try:
            to_pickle(self.formatted_results, filename)
        except:
            print("WARNING: COULD NOT SAVE")

    def test_load(self):
        obj = pickle.load(open(self.filename, "rb"))
        total_size, size_dict = size_of(obj)
        print("Total size (bytes): {}".format(total_size))
        return obj, total_size, size_dict