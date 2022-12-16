
import os
import glob

def mk_save_path(method_dir, dataset_dir, ext=None):
    if not os.path.exists(method_dir):
        os.mkdir(method_dir)

    method_dataset_dir = os.path.join(method_dir, dataset_dir)
    if not os.path.exists(method_dataset_dir):
        os.mkdir(method_dataset_dir)
    else:
        existing_dirs = glob.glob(method_dataset_dir + "*")
        method_dataset_dir += "_{}".format(len(existing_dirs))
        os.mkdir(method_dataset_dir)
    
    if not ext:
        use_ext = ""
    else:
        use_ext = ext
        
    result_path = os.path.join(method_dataset_dir, "benchmark_results{}".format(use_ext))
    
    if not ext:
        os.mkdir(result_path)
    
    return result_path