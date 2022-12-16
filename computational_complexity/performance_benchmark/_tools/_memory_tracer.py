
# -- import packages: --------------------------------------------------------------------
import psutil
import tracemalloc

# -- import local dependencies: ----------------------------------------------------------


# -- supporting function(s): -------------------------------------------------------------
def format_allocated_memory(scale=(1024 * 1024)):
    
    (curr, peak) = tracemalloc.get_traced_memory()
    return {"curr": curr/1e6, "peak": peak/1e6}


# -- MemoryUsage class: ------------------------------------------------------------------
class MemoryUsage(object):
    def __init__(self):
        self._per_mb = (1024 * 1024)
        
        
        self._alloc_mem = format_allocated_memory(self._per_mb)
        self._res_mem = psutil.Process().memory_info().rss  / self._per_mb
        
        tracemalloc.clear_traces()
        tracemalloc.stop()
        
    @property
    def resident(self):
        return self._res_mem
    
    @property
    def currently_allocated(self):
        return self._alloc_mem['curr']
    
    @property
    def peak_allocated(self):
        return self._alloc_mem['peak']
    
    @property
    def components(self):
        return {
            "resident": self._res_mem,
            "allocated_current": self._alloc_mem['curr'],
            "allocated_peak": self._alloc_mem['peak'],
               }
    
    def __repr__(self):
        
        identity = "Memory occupance (Mb):\n"
        identity += "-" * int(len(identity)-1) + "\n"
        for prop in self.__dir__():
            if (not prop.startswith("_")) and (prop != "components"):
                identity += "{}: {}\n".format(prop, getattr(self, prop))
        
        return identity
    
    
# -- memory tracer decorator: ------------------------------------------------------------
def memory_tracer(func):
    def wrapper(*args, **kwargs):
                
        tracemalloc.start()
        result = func(*args, **kwargs) 
        mu = MemoryUsage() 

        if not result is None:
            return result, mu
        return mu
        
    return wrapper