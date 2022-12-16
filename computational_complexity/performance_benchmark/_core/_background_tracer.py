
import os
import glob
import psutil
from datetime import datetime
import threading
import time

class BackgroundTracer:
    count = 0
    def __init__(self, path, quiet=False):
        
        self.path = path
        self.f = open(self.path, 'w')
        keys = ["datetime"] + list(psutil.virtual_memory()._asdict().keys()) + ["\n"]
        
        if not quiet:
            print("Writing memory trace to: {}".format(self.path))
        
        self.f.write(",".join(keys))
        self.f.close()
            
    def update(self):
        
        self.f = open(self.path, 'a')
        dt = datetime.now().isoformat()
        _vals = list(psutil.virtual_memory()._asdict().values())
        vals = [dt] + [str(v) for v in _vals] + ["\n"]
        self.f.write(",".join(vals))
        self.f.close()
        self.count += 1
        

# -- func: ---------------------------
def _memory_trace():
    
    path = "simba.background_memory_trace.txt"
        
    if os.path.exists(path):
        pname = path[:-4]
        count = len(glob.glob(pname + "*"))
        path  = "{}.{}.txt".format(pname, count)
    
    bt = BackgroundTracer(path)
    
    for _ in iter(int, 1):
        time.sleep(1)
        bt.update()
        
    return path
        
def _print_memory_tracing_thread(thread):
    
    print("\nMemory tracing background thread:")
    print("---------------------------------")
    print("Identity: {}".format(thread.ident))
    print("Name: {}".format(thread.name))
    print("Native ID: {}".format(thread.native_id))
        
def memory_trace():
    
    t = threading.Thread(name='child procs', target=_memory_trace)
    t.start()
    
    _print_memory_tracing_thread(t)    
        