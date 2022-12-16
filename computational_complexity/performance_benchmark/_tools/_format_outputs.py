

import licorice_font
import pickle

def to_pickle(obj, filename, protocol=pickle.HIGHEST_PROTOCOL, verbose=True):
    pickle.dump(obj, open(filename, "wb"), protocol=protocol)
    if verbose:
        msg = licorice_font.font_format("Successfully saved to", ["BOLD"])
        print("{}: {}".format(msg, filename))


def rm_adata_simba(runner):
    delattr(runner, "_adata")
    delattr(runner, "adata_CG")

def format_outputs_simba(outs, filename):

    runner, resource_df = outs[0][0]
    total_m, total_t = outs[0][1], outs[1]
    
    rm_adata_simba(runner)

    formatted = {
        "m": total_m,
        "t": total_t,
        "runner": runner,
        "resource_df": resource_df,
    }
    
    try:
        to_pickle(formatted, filename)
    except:
        print("WARNING: COULD NOT SAVE")
        
    return formatted


def rm_adata(runner):
    del runner._adata
    runner.outs = list(runner.outs)
    del runner.outs[0]
        
    return runner

def format_outputs_scanpy(outs, filename, include_runner=False):
    """ """
    runner, resource_df = outs
    runner = rm_adata(runner)
    
    formatted = {
        "m": runner.m,
        "t": runner.t,
        "runner": runner,
        "resource_df": resource_df,
    }
    to_pickle(formatted, filename)

    return formatted