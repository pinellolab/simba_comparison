
__module_name__ = "_traced_subroutine.py"
__doc__ = """Base class for the method-tracing sub-routine."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["mvinyard@broadinstitute.org"])


# -- import packages: --------------------------------------------------------------------
from abc import ABC, abstractmethod


# -- import local dependencies: ----------------------------------------------------------
from .._tools import time_tracer, memory_tracer, format_resource_dict, format_resource_df


# -- Base class: -------------------------------------------------------------------------
class TracedSubRoutine(ABC):
    __name__ = ""
    def __parse__(self, kwargs, ignore=["self"]):
        """Parse keyword arguments"""
        self._kwargs = {}
        for k, v in kwargs.items():
            if not k in ignore:
                self._kwargs[k] = v

    @time_tracer
    @memory_tracer
    def __tracer__(self, kwargs):
        """Trace forward"""
        return self.forward(**kwargs)

    def __call__(self, **kwargs):
        """Call the forward method, which is wrapped by the tracer"""
        self.__parse__(locals())
        result = self.__tracer__(**self._kwargs)
        self.outs = result[0][0]
        self.m, self.t = result[0][1], result[1]
        self._ResourceDict = format_resource_dict(self, self.__name__)
        self.resource_df = format_resource_df(self._ResourceDict)

    @abstractmethod
    def forward(self):
        """
        This method should be overwritten / customized. This is where the
        an analysis method sub-routine is run; this function is traced.
        """
        pass