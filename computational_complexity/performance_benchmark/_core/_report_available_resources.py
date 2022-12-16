import psutil
import os

class Memory:
    def __init__(self, use_mb=False, ignore=['percent', 'index', 'count']):
        
        self.svmem = psutil.virtual_memory()
        
        for attr in self._attributes():
            if not attr in ignore:
                mem_attr = getattr(self.svmem, attr)
                if not use_mb:
                    mem_attr = "{} Gb".format(self._to_Gb(mem_attr))
                else:
                    mem_attr = "{} Mb".format(self._to_Mb(mem_attr))
                setattr(self, attr, mem_attr)
        
    def _attributes(self):
        return [i for i in self.svmem.__dir__() if not i.startswith("_")]
    
    def _to_Gb(self, val):
        return val >> 30

    def _to_Mb(self, val):
        return val >> 20
    


def _available_memory(use_mb=False):
    """report the available memory in Mb or Gb"""
    total = Memory(use_mb=use_mb).total
    return "Available memory: {}".format(total)

def _available_cpus():
    return "Available vCPUs: {}".format(os.cpu_count())

def report_available_resources():
    print(_available_memory())
    print(_available_cpus())