
import simba as si


from ...._core import TracedSubRoutine
from .... import _tools as tl

memory_tracer = tl.memory_tracer
time_tracer = tl.time_tracer


class SimbaPreprocess(TracedSubRoutine):
    
    __name__ = "simba_pp"

    _functions = [
        "pp.filter_genes",
        "pp.cal_qc_rna",
        "pp.filter_cells_rna",
        "pp.normalize",
        "pp.log_transform",
        "pp.select_variable_genes",
        "tl.discretize",
    ]

    def forward(
        self,
        adata,
        min_n_cells=3,
        min_n_genes=200,
        quiet=True,
        method="lib_size",
        n_top_genes=2000,
        n_bins=5,
    ):

        kwargs = locals()

        adata_CG = adata.copy()
        ResourceDict = {"m": {}, "t": {}}

        for mod_func in self._functions:
            mod, func = mod_func.split(".")

            @memory_tracer
            @time_tracer
            def f(run, adata_CG, kwargs):
                kw = tl.extract_func_kwargs(run, kwargs)
                kw["adata"] = adata_CG
                return run(**kw)

            t, m = f(getattr(getattr(si, mod), func), adata_CG, kwargs)
            ResourceDict["m"][mod_func] = m.components
            ResourceDict["t"][mod_func] = t

        return adata_CG, tl.format_resource_df(ResourceDict)