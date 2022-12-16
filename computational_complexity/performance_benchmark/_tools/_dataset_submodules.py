
from abc import ABC

class MethodModules(ABC):
    @property
    def properties(self):
        self._properties = [attr for attr in self.__dir__() if not attr.startswith("_")]
        self._properties.remove("properties")
        return self._properties


class SIMBAModules(MethodModules):
    @property
    def pp_io(self):
        self._pp = [
            "pp.filter_genes",
            "pp.cal_qc_rna",
            "pp.filter_cells_rna",
            "pp.normalize",
            "pp.log_transform",
            "pp.select_variable_genes",
            "tl.discretize",
            "modify_params",
            "edge_write_time"
        ]
        return self._pp

    @property
    def graph_gen(self):
        self._gen_graph = ["graph_gen_time"]
        return self._gen_graph

    @property
    def train(self):
        self._train = ["train"]
        return self._train

    @property
    def embedding(self):
        self._emb = [
            "read_embed",
            "cell_umaps",
            "joint_embed",
            "annot_entities",
            "joint_umap",
        ]
        return self._emb

    @property
    def gene_discovery(self):
        self._gene_disc = ["compare_entities"]
        return self._gene_disc


class ScanpyModules(MethodModules):
    @property
    def pp(self):
        self._pp = [
            "pp.filter_cells",
            "pp.filter_genes",
            "annot_mt",
            "pp.calculate_qc_metrics",
            "do_qc_filter",
            "pp.normalize_total",
            "pp.log1p",
            "pp.highly_variable_genes",
            "do_hv_filter",
            "pp.regress_out",
            "pp.scale",
        ]

        return self._pp

    @property
    def embedding(self):
        self._emb = [
            "tl.pca",
            "pp.neighbors",
            "tl.umap",
        ]
        return self._emb

    @property
    def gene_discovery(self):
        self._gene_disc = ["tl.leiden", "tl.rank_genes_groups"]
        return self._gene_disc


class SeuratModules(MethodModules):
    @property
    def pp(self):
        self._pp = [
            "CreateSeuratObject",
            "PercentageFeatureSet",
            "NormalizeData",
            "FindVariableFeatures",
            "ScaleData",
        ]

        return self._pp

    @property
    def embedding(self):
        self._emb = ["RunPCA", "FindNeighbors", "RunUMAP"]
        return self._emb

    @property
    def gene_discovery(self):
        self._gene_disc = ["FindClusters", "FindAllMarkers"]
        return self._gene_disc
    
    
class MethodSubModules:
    def __init__(self):
        self.SIMBA = SIMBAModules()
        self.Scanpy = ScanpyModules()
        self.Seurat = SeuratModules()