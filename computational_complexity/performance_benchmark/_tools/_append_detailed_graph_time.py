
import pandas as pd

def _extract_detailed_graph_time(results):
    
    gen_graph_resource_df = results['runner']._gen_graph.resource_df
    outs_dict = results['runner']._gen_graph.outs
    outs_dict["graph_gen_time"] = gen_graph_resource_df.loc['gen_graph_total', 't'] - outs_dict['edge_write_time']
    tmp_df = pd.DataFrame.from_dict(outs_dict, orient="index").rename({0:"t"}, axis=1)
    
    return pd.concat([gen_graph_resource_df, tmp_df])

def append_detailed_graph_time(results):
    return pd.concat([results['resource_df'], _extract_detailed_graph_time(results)])
