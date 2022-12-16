

import os
import simba as si


from ...._core import TracedSubRoutine


class ModifyParameters(TracedSubRoutine):
    __name__ = "modify_params"
    def forward(self, dataset):
        dict_config = si.settings.pbg_params.copy()
        dict_config["workers"] = os.cpu_count()
        return dict_config