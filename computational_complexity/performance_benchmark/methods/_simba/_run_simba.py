

# -- import packages: ---------------
from tqdm.notebook import tqdm
from datetime import datetime
import pandas as pd
import simba as si
import glob
import os

# -- import local dependencies: ---------------
from . import _simba_funcs as funcs

from ... import _tools as tl
    
time_tracer = tl.time_tracer
memory_tracer = tl.memory_tracer


from ._format_simba_results import format_simba_results


# -- supporting functions: ------------------------------
def _n_cells(adata):
    return adata.shape[0]

def format_graph_time(runner):
#     run_i[0][0]
    resource_df = runner[1]
    gen_graph = runner[0]._gen_graph

    graph_total_t = gen_graph.resource_df["t"]
    edge_write_t  = gen_graph.outs["edge_write_time"]

    g_df = pd.concat(
        [
            gen_graph.resource_df,
            pd.Series(graph_total_t - edge_write_t).to_frame(),
        ],
    )
    g_df.index = ["write_graph", "gen_graph"]
    
    print("\n\nresource_df")
    print(resource_df)
    print("\n\ng_df")
    print(g_df)

    return pd.concat([resource_df, g_df])


def _format_resources(resources):
    
    for key in ["train", "gen_graph_total"]: # "gen_graph_detailed"
        resources[key].index = [key]

    return pd.concat(
        [
            resources["preproc_detailed"],
            resources["mod_params"],
            resources["train"],
            resources["postproc"],
        ]
    )    
#     format_graph_time(runner, resource_df)
    

    
# -- Method operating class: ------------------------------
class RunSIMBA:
    def __init__(self, adata):
        self._adata = adata

    def preprocess(self):
        self._simba_pp = funcs.SimbaPreprocess()
        self._simba_pp(adata=self._adata)
        self.adata_CG = self._simba_pp.outs[0]

    def gen_graph(self, dataset, use_hv=True):

        self.dataset = dataset

        self._gen_graph = funcs.GenerateGraph()
        self._gen_graph(adata_CG=self.adata_CG, use_hv=use_hv, dataset=dataset)

    def modify_params(self):
        self._modify_params = funcs.ModifyParameters()
        self._modify_params(dataset=self.dataset)

    def train(self):

        self._simba_train = funcs.SimbaTrain()
        self._simba_train(dataset=self.dataset, dict_config=self._modify_params.outs)
        self.params = si.settings.pbg_params

    def post_process(
        self,
    ):

        self.simba_post_process = funcs.SimbaPostProcess(params=self.params)
        self.simba_post_process()
        
    def collect_resources(self):
        
        self.resources = {}
        
        self.resources['preproc_detailed'] = self._simba_pp.outs[1]
        self.resources['preproc_total'] = self._simba_pp.resource_df
        self.resources['gen_graph_total'] = self._gen_graph.resource_df
#         self.resources['gen_graph_detailed'] = pd.DataFrame.from_dict(
#             self._gen_graph.outs[1], orient="index"
#         ).T
        self.resources['train'] = self._simba_train.resource_df
        self.resources['mod_params'] = self._modify_params.resource_df
        self.resources['postproc'] = self.simba_post_process.post_process_df
        
        self.resource_df = _format_resources(self.resources)
        
        return self.resource_df
    
# -- trace execution of the operating class: ------------------------------
def now():
    return datetime.now().isoformat()
    


@time_tracer
@memory_tracer
def _run_simba_traced(adata, dataset):
    
    t = {'t0': now()}

    run_simba = RunSIMBA(adata)
    t['init'] = now()
    
    run_simba.preprocess()
    t['pp'] = now()
    
    run_simba.gen_graph(dataset=dataset)
    t['gen_graph'] = now()
    
    run_simba.modify_params()
    t['modify_params'] = now()
    
    run_simba.train()
    t['train'] = now()
    
    run_simba.post_process()
    t['post_process'] = now()
    
    resource_df = run_simba.collect_resources()
    t['tf'] = now()
        
    return run_simba, resource_df, t

# -- API-facing exe and collection of results: ------------------------------
def run_simba(adata,
              n_iters=5,
              method_dir="./simba_benchmark_results",
              dataset="unspecified_dataset",
              parent_path=os.getcwd(),
              graph_dir="./result_simba/pbg",
              n_cells=None,
              debug=False,
             ):
    
    outs = {}
    
    if not n_cells:
        n_cells = _n_cells(adata)
    
    for i in tqdm(range(n_iters)):    
        filename = tl.mk_save_path(method_dir, dataset, ext=".pkl")
        tl.clean_graph_dir(dataset=dataset, parent_path=parent_path, graph_dir=graph_dir)
        run_i = _run_simba_traced(adata, dataset)
     
        outs[i] = format_simba_results(
            run=run_i,
            save_filename=filename,
            dataset=dataset,
            n_cells=n_cells,
            run_iter=i,
            include_runner=False,
            return_saved_outs=True,
        )
        if debug:
            return run_i, outs[i]

    return outs
 