
# -- import packages: --------------------------------------------------------------------
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import vinplots
import os


# -- import local dependencies: ----------------------------------------------------------
from ._labeled_color import LabeledColor
from ._method_functions import MethodFunctions
from ._result_classes import SIMBAResults, ScanpyResults, SeuratResults
from ._fetch_properties import fetch_properties


# -- supporting functions: ---------------------------------------------------------------
def function_group_resources(df, method_class):
    grouped = {}
    for prop in fetch_properties(method_class):
        grouped[prop] = df.loc[getattr(method_class, prop)].sum()

    return grouped


def parse_results(results_obj, method_funcs, seurat=False):

    collect = Collect(results_obj())

    if not seurat:
        time_df = collect.time().mean(1)
        mem_df  = collect.memory().mean(1)
    else:
        time_df = collect.time(seurat=True).mean(1)
        mem_df  = collect.memory(seurat=True).mean(1)

    func_time = function_group_resources(time_df, method_funcs)
    func_mem  = function_group_resources(mem_df,  method_funcs)
    
    return {"t": func_time, "m": func_mem}

        
def cumulative_barplot(ax, df, labcol, log=True, legend=True, legend_loc=(2.5,0)):

    ax.set_title("Time", fontsize=10)

    for n, (i, row) in enumerate(df.iterrows()):
        if n > 0:
            bottom = df.iloc[int(n - 1)]
        else:
            bottom = 0
        ax.bar(
            x=range(3),
            height=row,
            bottom=bottom,
            label=labcol.loc[i]["labels"],
            color=labcol.loc[i]["colors"],
            edgecolor="k",
            log=log,
        )

    ax.set_xticks([0, 1, 2], ["SIMBA", "Scanpy", "Seurat"], fontsize=8)
    ax.set_ylabel("s", rotation=0, fontsize=8)

    ax.tick_params(axis="both", which="major", labelsize=6)
    if legend:
        ax.legend(edgecolor="w", loc=legend_loc)


def memory_plot(
    ax,
    df,
    labcol,
    space=1.5,
    mem_keys=["pp", "generate_graph", "embedding", "train", "gene_discovery"],
    log=True,
    legend=False,
):

    ax.set_title("Memory", fontsize=10)

    plot_mem = df.loc[mem_keys]
    count = 0
    xr_ = np.array([])
    for n, col in enumerate(plot_mem.columns):
        n_rows = plot_mem[col].shape[0]
        ax.bar(
            x=np.arange(count, count + n_rows),
            height=plot_mem[col],
            color=labcol.loc[mem_keys]["colors"],
            label=col,
            edgecolor="k",
            lw=0.5,
            log=log,
        )
        count += n_rows + space

    ax.set_xticks([2, 8.5, 15], ["SIMBA", "Scanpy", "Seurat"], fontsize=8)
    ax.tick_params(axis="both", which="major", labelsize=6)
    ax.set_ylabel("Mb", rotation=0, fontsize=8)

    if legend:
        ax.legend(edgecolor="w")

        
def process_results_for_plot(compare, tm_key, idx=["SIMBA", "Scanpy", "Seurat"]):

    results = [compare.Results[key][tm_key] for key in idx]
    df = pd.DataFrame(results, index=idx).T.fillna(0).fillna(0)
    if tm_key == "t":
        return df.cumsum()
    return df


class Collect:
    def __init__(self, ResultsDict):
        self.ResultsDict = ResultsDict

    def time(self, seurat=False):

        Time = {}
        for key in self.ResultsDict.keys():
            if not seurat:
                Time[key] = self.ResultsDict[key]["resource_df"]["t"]
            else:
                Time[key] = self.ResultsDict[key].time["time_dict"]

        return pd.DataFrame(Time)

    def memory(self, seurat=False):

        Memory = {}
        for key in self.ResultsDict.keys():
            if not seurat:
                Memory[key] = self.ResultsDict[key]["resource_df"]["allocated_peak"]
            else:
                Memory[key] = self.ResultsDict[key].compute_mb_used()

        return pd.DataFrame(Memory)

    
class CompareMethods:
    
    method_funcs = MethodFunctions()
    
    _method_paths = {
        "SIMBA": "SIMBA/*.pkl",
        "Scanpy": "Scanpy/*.pkl",
        "Seurat": "Seurat/*/",
    }

    def __init__(self, path):
        self._base_path = path
        for method, m_path in self._method_paths.items():
            method_path = os.path.join(path, m_path)
            setattr(self, "_{}_path".format(method), method_path)

        self.Results = {}
        self.Results["SIMBA"] = parse_results(
            SIMBAResults(self._SIMBA_path), self.method_funcs.SIMBA
        )
        self.Results["Scanpy"] = parse_results(
            ScanpyResults(self._Scanpy_path), self.method_funcs.Scanpy
        )
        try:
            self.Results["Seurat"] = parse_results(
                SeuratResults(self._Seurat_path), self.method_funcs.Seurat, seurat=True
            )
        except:
            self.Results["Seurat"] = {"t": {"pp":0, "embedding":0, "gene_discovery" :0}, 
                                      "m": {"pp":0, "embedding":0, "gene_discovery" :0},
                                     }
        self._SIMBA_RESULTS = SIMBAResults(self._SIMBA_path)


# -- Augment SIMBA memory usage from trace file: -----------------------------------------
class SIMBAMemoryAugmentation:
    def __init__(self, compare):
        
        self.simba_outs = compare._SIMBA_RESULTS()
    
    def _time_from_iso(self, t):
        s = datetime.datetime.fromisoformat(t).second / (60*60)
        m = datetime.datetime.fromisoformat(t).minute / 60
        h = datetime.datetime.fromisoformat(t).hour
        
        return h + m + s
    
    def parse(self):
        
        self.MemDict = {}
        
        for k, v in self.simba_outs.items():
            self.MemDict[k] = {}
            self.MemDict[k]['landmarks'] = v['landmark_datetimes']
            self.MemDict[k]['memory_df'] = v['memory_trace_df']

    def _filter_time_df(self, time_df, landmark_dict, t0_key, tf_key):

        t0_ = landmark_dict[t0_key]
        tf_ = landmark_dict[tf_key]

        return time_df.loc[time_df['datetime'] > t0_].loc[time_df['datetime'] < tf_]
            
    def resolve_time(self, landmarks, time_df):
    
        keys = list(landmarks.keys())
        min_active = (time_df['active'].min() / 1e6)

        ResolvedTime = {}
        for i, key in enumerate(keys):
            if not i == 0:
                _df = self._filter_time_df(time_df, landmarks, keys[i-1], key) 
                if _df.shape[0] > 0:
                    # some steps are so fast, the memory is not even recorded
                    ResolvedTime[key] = _df

        MaxMemAugmentation = {}
        for key, _df in ResolvedTime.items():
            max_mem =  _df['used'].max()/1e6 - min_active
#             print("{:<15}{:.2f} Mb".format(key,max_mem))
            MaxMemAugmentation[key] = max_mem
        return MaxMemAugmentation
    
    def update(self, plot_m):
        
        ALL = {}
        self.parse()
        for i in self.MemDict.keys():
            landmarks = self.MemDict[i]['landmarks']
            time_df = self.MemDict[i]['memory_df']
            ALL[i] = self.resolve_time(landmarks, time_df)
            
            
        augment_df = pd.DataFrame(ALL).loc[['gen_graph', 'train']]
        cols = augment_df.columns
        if len(cols) > 3:
            cols = cols[1:]
        augmented = augment_df[cols].mean(1).to_dict()
        
        
        plot_m.loc['train']['SIMBA'] = augmented['train']
        plot_m.loc['generate_graph']['SIMBA'] = augmented['gen_graph']
        
        return plot_m
    
# -- Main module class: ------------------------------------------------------------------
class PlotTimeMemory:
    def __init__(self,
                 results_path,
                 memory_trace_path=None,
                 wspace=0.2,
                 save=False,
                 savename=".svg",
                 legend_loc=(1,0),
                 **kwargs,
                ):

        # -- gather results: -------------------------------------------------------------
        self._compare = CompareMethods(results_path)
        self._results_path = results_path
        SIMBA_mem_augment = SIMBAMemoryAugmentation(self._compare)

        self.plot_m = process_results_for_plot(self._compare, tm_key="m")
        self.plot_m = SIMBA_mem_augment.update(self.plot_m)
        
        self.plot_t = process_results_for_plot(self._compare, tm_key="t")
        
        
        # -- plot: -----------------------------------------------------------------------
        self._labeled_colors = LabeledColor().frame()
        self.fig, self._axes = vinplots.quick_plot(nplots=2, ncols=2, wspace=wspace, **kwargs)
        self.t_ax, self.m_ax = self._axes[0], self._axes[1]
        
        cumulative_barplot(self.t_ax, self.plot_t, labcol=self._labeled_colors, legend=True, legend_loc=legend_loc)
        memory_plot(self.m_ax, self.plot_m, labcol=self._labeled_colors, legend=False)
        
        if save:
            plt.savefig(savename)
