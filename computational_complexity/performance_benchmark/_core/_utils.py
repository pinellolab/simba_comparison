
from datetime import datetime
import pandas as pd
import glob
import os


# -- funcs: ---------------------------------------------
def time_from_iso(t):
    s = datetime.fromisoformat(t).second / (60*60)
    m = datetime.fromisoformat(t).minute / 60
    h = datetime.fromisoformat(t).hour
    return h + m + s





def read_latest_memory_trace(path="simba.background_memory_trace*.txt"):
    
    memory_trace_files = {}
    
    for p in glob.glob(path):
        split_p = os.path.basename(p).split(".")
        if not len(split_p) == 4:
            version = 0
        else:
            version = split_p[-2]
        memory_trace_files[int(version)] = p
        
    mem_paths = pd.DataFrame.from_dict(memory_trace_files, orient="index").reset_index()
    mem_paths.columns = ["iter", "path"]
    mem_paths['iter'] = mem_paths['iter'].astype(int)
    
    mem_path_df = mem_paths.sort_values('iter').reset_index(drop=True)
    
    use_path = mem_path_df.iloc[mem_path_df['iter'].iloc[-1]]['path']
    mem_df = pd.read_csv(use_path, usecols=range(12))
    
    return mem_df, use_path