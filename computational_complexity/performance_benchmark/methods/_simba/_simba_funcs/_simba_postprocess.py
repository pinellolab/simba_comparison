
# -- : -----------------------------------------------------------------------------
import simba as si
import pandas as pd
import os


# -- : -----------------------------------------------------------------------------
from ...._core import TracedSubRoutine


# -- : -----------------------------------------------------------------------------
class ReadEmbeddings(TracedSubRoutine):
    def forward(self): # params
        """
        C: embeddings for cells
        G: embeddings for genes
        """
        dict_adata = si.read_embedding() # path_emb=params["checkpoint_path"], 
#                                  path_entity=params["entity_path"],
#                                 )
        return dict_adata

class CellUMAP(TracedSubRoutine):
    """returns: None, modifies adata_C"""

    def forward(self, adata_C, n_neighbors=15, n_components=2):
        si.tl.umap(adata_C, n_neighbors=n_neighbors, n_components=n_components)
        return "NULL"

class JointEmbed(TracedSubRoutine):
    def forward(self, adata_C, adata_G):
        """returns: adata_all"""
        return si.tl.embed(adata_ref=adata_C, list_adata_query=[adata_G])


class AnnotateEntities(TracedSubRoutine):
    def forward(self, adata_all, adata_G, entity="gene"):
        """returns adata_all"""
        adata_all.obs["entity_anno"] = ""
        adata_all.obs.loc[adata_G.obs_names.astype(str), "entity_anno"] = entity
        return adata_all


class JointUMAP(TracedSubRoutine):
    """returns: None, modifies adata_all"""

    def forward(self, adata_all, n_neighbors=15, n_components=2):
        si.tl.umap(adata_all, n_neighbors=n_neighbors, n_components=n_components)
        return "NULL"

class CompareEntities(TracedSubRoutine):
    def forward(self, adata_C, adata_G):
        """returns adata_cmp"""
        return si.tl.compare_entities(adata_ref=adata_C, adata_query=adata_G)


# -- : -----------------------------------------------------------------------------
class SimbaPostProcess:
    def __init__(self, params):

        self.params = params
        
        self.read_embed = ReadEmbeddings()
        self.cell_umaps = CellUMAP()
        self.joint_embed = JointEmbed()
        self.annot_entities = AnnotateEntities()
        self.joint_umap = JointUMAP()
        self.compare_entities = CompareEntities()

    def format_post_process_resources(self):
    
        funcs = [
            "read_embed",
            "cell_umaps",
            "joint_embed",
            "annot_entities",
            "joint_umap",
            "compare_entities",
        ]

        self.post_process_df = pd.concat(
            [getattr(self, func).resource_df for func in funcs]
        )
        self.post_process_df.index = funcs
        
    def __call__(self, n_neighbors=15, n_components=2, entity="gene"):
        
        self.read_embed()
        self.data_dict = self.read_embed.outs
        adata_C, adata_G = self.data_dict["C"], self.data_dict["G"]
        self.cell_umaps(adata_C=adata_C, n_neighbors=n_neighbors, n_components=n_components)
        self.joint_embed(adata_C=adata_C, adata_G=adata_G)
        adata_all = self.joint_embed.outs
        self.annot_entities(adata_all=adata_all, adata_G=adata_G, entity=entity)
        adata_all = self.annot_entities.outs
        self.joint_umap(adata_all=adata_all, n_neighbors=n_neighbors, n_components=n_components)
        self.compare_entities(adata_C=adata_C, adata_G=adata_G)
        
        self.adata_C = adata_C
        self.adata_G = adata_G
        self.adata_all = adata_all
        self.adata_cmp = self.compare_entities.outs

        self.format_post_process_resources()