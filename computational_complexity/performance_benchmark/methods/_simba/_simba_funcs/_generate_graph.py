
from ...._core import TracedSubRoutine

import simba as si

class GenerateGraph(TracedSubRoutine):
    def forward(self, adata_CG, use_hv, dataset):
        _, edge_write_time = si.tl.gen_graph(
            list_CG=[adata_CG], copy=False, use_highly_variable=use_hv, dirname=dataset
        )
        return edge_write_time