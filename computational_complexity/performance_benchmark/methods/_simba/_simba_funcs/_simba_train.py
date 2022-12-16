

# -- training: ----------------------------------------------------------------------------------------------------
import simba as si


# -- training: ----------------------------------------------------------------------------------------------------
from ...._core import TracedSubRoutine
from ...._tools import time_tracer, memory_tracer


# -- training: ----------------------------------------------------------------------------------------------------
class SimbaTrain(TracedSubRoutine):
    def __init__(self):
        pass

    @time_tracer
    @memory_tracer
    def forward(self, dataset, dict_config):
        
        return si.tl.pbg_train(
            dirname=dataset,
            pbg_params=dict_config,
            auto_wd=True,
            save_wd=True,
            output="model",
        )