import scanpy as sc

from ..._core import TracedSubRoutine

from ... import _tools as tl
    
time_tracer = tl.time_tracer
memory_tracer = tl.memory_tracer


# -- supporting functions: ------------------------------
def do_qc_filter(adata, max_pct_mt, max_n_genes_by_counts):
    adata = adata[adata.obs.n_genes_by_counts < max_n_genes_by_counts, :]
    adata = adata[adata.obs.pct_counts_mt < max_pct_mt, :]
    return adata.copy()


def do_hv_filter(adata):
    adata.raw = adata
    return adata[:, adata.var.highly_variable].copy()


def annot_mt(adata):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")


# -- key class: ---------------------------------------------
class TracedScanpy(TracedSubRoutine):

    __name__ = "scanpy"

    _functions = [
        "pp.filter_cells",
        "pp.filter_genes",
        "annot_mt",
        "pp.calculate_qc_metrics",
        "do_qc_filter",
        "pp.normalize_total",
        "pp.log1p",
        "pp.highly_variable_genes",
        "do_hv_filter",
        "pp.regress_out",
        "pp.scale",
        "tl.pca",
        "pp.neighbors",
        "tl.umap",
        "tl.leiden",
        "tl.rank_genes_groups",
    ]

    def forward(
        self,
        adata,
        dataset=None,
        do_plots=False,
        inplace=True,
        percent_top=None,
        log1p=False,
        min_genes=200,
        min_cells=3,
        max_n_genes_by_counts=2000,
        max_pct_mt=5,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_neighbors=10,
        n_pcs=40,
        qc_vars=["mt"],
        target_sum=1e4,
        keys=["total_counts", "pct_counts_mt"],
        max_value=10,
        svd_solver="arpack",
        groupby="leiden",
        method="t-test",
        nn_method="umap",
    ):

        kwargs = locals()

        self._adata = adata
        ResourceDict = {"m": {}, "t": {}}

        for mod_func in self._functions:
            if "." in mod_func:
                mod, func = mod_func.split(".")
                func_to_run = getattr(getattr(sc, mod), func)
            else:
                func_to_run = globals()[mod_func]

            @memory_tracer
            @time_tracer
            def f(run, adata, kwargs):
                args = []
                kw = tl.extract_func_kwargs(run, kwargs)
                if run.__str__().split(" ")[1] in [
                    "filter_cells",
                    "filter_genes",
                    "pca",
                ]:
                    kw["data"] = adata

                elif run.__str__().split(" ")[1] in ["log1p", "scale"]:
                    args.append(adata)
                else:
                    kw["adata"] = adata

                if run.__str__().split(" ")[1] in ["umap", "neighbors"]:
                    kw["method"] = nn_method

                return run(*args, **kw)

            if func_to_run.__str__().split(" ")[1] in ["do_qc_filter", "do_hv_filter"]:
                (self._adata, t), m = f(func_to_run, self._adata, kwargs)
            else:
                t, m = f(func_to_run, self._adata, kwargs)
            ResourceDict["m"][mod_func] = m.components
            ResourceDict["t"][mod_func] = t
            
        return self._adata, tl.format_resource_df(ResourceDict)