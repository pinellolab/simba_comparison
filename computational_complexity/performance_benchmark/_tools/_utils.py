import inspect
import pandas as pd


def func_params(func):
    return list(inspect.signature(func).parameters.keys())


def extract_func_kwargs(func, kwargs):
    func_kwargs = {}
    params = func_params(func)
    for k, v in kwargs.items():
        if k in params:
            func_kwargs[k] = v
    return func_kwargs


def _parse_kwargs(kwargs, ignore=["adata"]):

    _kwargs = {}

    for k, v in kwargs.items():
        if not k in ignore:
            _kwargs[k] = v
    return _kwargs


def format_resource_df(ResourceDict):

    t_df = pd.DataFrame.from_dict(ResourceDict["t"], orient="index", columns=["t"])
    m_df = pd.DataFrame(ResourceDict["m"]).T

    return pd.concat([m_df, t_df], axis=1)


def format_resource_dict(routine, name="routine"):
    """format the t/m resource routine"""
    return {
        "t": {name: routine.t},
        "m": {name: routine.m.components},
    }