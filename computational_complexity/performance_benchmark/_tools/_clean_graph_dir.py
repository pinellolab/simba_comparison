

import os

def clean_graph_dir(dataset, parent_path=os.getcwd(), graph_dir="./result_simba/pbg"):
    """dataset must match the dataset name given to simba - this way, we remove the exactly correct directory.
    To do this, we pass the filename created by tl.mk_save_path and get the immediate parent dir.
    """
        
    path = os.path.join(parent_path, graph_dir, dataset)
    if os.path.exists(path):
        os.system("rm -r {}".format(path))
        print(" - Removed graph dir: {}".format(path))