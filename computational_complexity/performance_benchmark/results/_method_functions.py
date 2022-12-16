from abc import ABC


class CategoricalFunctions(ABC):
    def __init__(self):

        self._pp = []
        self._embedding = []
        self._gene_discovery = []
        self._generate_graph = []
        self._write_graph = []
        self._train = []

    @property
    def pp(self):
        return self._pp

    @property
    def generate_graph(self):
        return self._generate_graph

    @property
    def write_graph(self):
        return self._write_graph

    @property
    def train(self):
        return self._train

    @property
    def embedding(self):
        return self._embedding

    @property
    def gene_discovery(self):
        return self._gene_discovery


class SIMBAFunctions(CategoricalFunctions):
    def __init__(self):

        self._pp = [
            "pp.filter_genes",
            "pp.cal_qc_rna",
            "pp.filter_cells_rna",
            "pp.normalize",
            "pp.log_transform",
            "pp.select_variable_genes",
            "tl.discretize",
            "modify_params",
            "read_embed",
            "annot_entities",
        ]
        self._generate_graph = ["gen_graph"]
        self._write_graph = ["write_graph"]
        self._train = ["train"]
        self._embedding = ["cell_umaps", "joint_embed", "joint_umap"]
        self._gene_discovery = ["compare_entities"]


class ScanpyFunctions(CategoricalFunctions):
    def __init__(self):

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
        self._embedding = [
            "tl.pca",
            "pp.neighbors",
            "tl.umap",
        ]
        self._gene_discovery = ["tl.leiden", "tl.rank_genes_groups"]


class SeuratFunctions(CategoricalFunctions):
    def __init__(self):
        self._pp = [
            "CreateSeuratObject",
            "PercentageFeatureSet",
            "NormalizeData",
            "FindVariableFeatures",
            "ScaleData",
        ]
        self._embedding = ["RunPCA", "FindNeighbors", "RunUMAP"]
        self._gene_discovery = ["FindClusters", "FindAllMarkers"]
        
        
class MethodFunctions:
    def __init__(self):
        pass

    @property
    def SIMBA(self):
        return SIMBAFunctions()
    
    @property
    def Scanpy(self):
        return ScanpyFunctions()
    
    @property
    def Seurat(self):
        return SeuratFunctions()