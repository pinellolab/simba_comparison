import time

def time_tracer(abstract_func):
    
    def timing_wrapper(*args, **kwargs):
        
        t0 = time.time()
        result = abstract_func(*args, **kwargs)
        tf = time.time()
        
        t_delta = tf - t0
        
        if not result is None:
            return result, t_delta
        return t_delta
    
    return timing_wrapper