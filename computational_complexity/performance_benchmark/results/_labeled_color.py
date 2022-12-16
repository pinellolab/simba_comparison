
import numpy as np
import pandas as pd


from ._fetch_properties import fetch_properties


class LabeledColor:
    def __init__(self):
        self.keys = [
            "pp",
            "generate_graph",
            "write_graph",
            "train",
            "embedding",
            "gene_discovery",
        ]
        self.colors = ["#eb5133", "#ff8133", "#ffc247", "#76cd65", "#88dae7", "#0d3082"]
        self.labels = [
            "Preprocessing",
            "Graph Generation",
            "Graph IO",
            "Training",
            "Embedding",
            "Gene Discovery",
        ]

    def frame(self):
        LC = {}
        for key in fetch_properties(self, ignore=["frame"]):
            LC[key] = getattr(self, key)

        return pd.DataFrame(LC).set_index("keys")