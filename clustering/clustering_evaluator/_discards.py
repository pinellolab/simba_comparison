# def evaluate(adata, savedir="./"):

#     clustering_outpath = os.path.join(savedir, adata.uns["method"] + "_clusters.tsv")
#     path_dir = "/".join(clustering_outpath.split("/")[:-1])
#     if os.path.exists(path_dir) == False:
#         os.mkdir(path_dir)

#     adata.obs[["louvain", "kmeans", "hc"]].to_csv(clustering_outpath, sep="\t")

#     metrics_outpath = os.path.join(savedir, adata.uns["method"] + "_metrics.tsv")
#     df_metrics.to_csv(metrics_outpath, sep="\t")

#     try:
#         adata.obs = adata.obs.reset_index()
#     except:
#         pass

#     if clustering == None:
#         celltypes = adata.obs.celltype.values.unique().astype(str)
#         for i in celltypes:
#             cells = adata.obs.loc[adata.obs.celltype == str(i)].index.astype(int)
#             x, y = adata.obsm['X_umap'][cells][:,0], adata.obsm['X_umap'][cells][:,1]
#             ax.scatter(x, y, s=4, c=adata.obs.coloring[cells], label=i)
#     elif clustering == "louvain":
#         clusters = adata.obs.louvain.unique()
#     elif clustering == "kmeans":
#         clusters = adata.obs.kmeans.unique()
#     elif clustering == "hc":
#         clusters = adata.obs.hc.unique()
#     else:
#         print("Please select a clustering type or empty the field to indicate None.")

#     if clustering != None:
#         for clust in np.sort(clusters):
#             cells = adata.obs.loc[adata.obs[clustering] == clust].index.astype(int)
#             x, y = adata.obsm['X_umap'][cells][:,0], adata.obsm['X_umap'][cells][:,1]
#             ax.scatter(x, y, s=4, c=extended_color_pallette[clust], label=clust)

#     ax.spines['left'].set_linewidth(3)
#     ax.spines['bottom'].set_linewidth(3)
#     ax.spines['right'].set_visible(False)
#     ax.spines['top'].set_visible(False)
#     ax.spines['left'].set_color('black')
#     ax.spines['bottom'].set_color('black')
#     ax.yaxis.set_ticks_position('left')
#     ax.xaxis.set_ticks_position('bottom')
#     ax.grid(b=None)
#     plt.legend(
#         markerscale=3,
#         edgecolor="w",
#         fontsize=14,
#         handletextpad=None,
#         bbox_to_anchor=(0.5, 0.0, 0.70, 1),
#     )
#     plt.tight_layout()

#     if save != False:
#         plt.savefig(save)

#     plt.legend(
#         markerscale=3,
#         edgecolor="w",
#         fontsize=14,
#         handletextpad=None,
#         bbox_to_anchor=(0.5, 0.0, 0.70, 1),
#     )
#     plt.show()


# """
# Eventually, use this color pallette:

# cbPalette <- c()

# Raw: "#999999", SIMBA: "#009E73"
# """

# class EmbeddingEvaluation:
#     def __init__(self):
#         pass

#     def clustering(self):
#         pass

#     def run_metrics(self):
#         pass

#     #adjusted rank index
#     ari_louvain = adjusted_rand_score(adata.obs['celltype'], adata.obs['louvain'])
#     ari_kmeans = adjusted_rand_score(adata.obs['celltype'], adata.obs['kmeans'])
#     ari_hc = adjusted_rand_score(adata.obs['celltype'], adata.obs['hc'])
#     #adjusted mutual information
#     ami_louvain = adjusted_mutual_info_score(adata.obs['celltype'], adata.obs['louvain'],average_method='arithmetic')
#     ami_kmeans = adjusted_mutual_info_score(adata.obs['celltype'], adata.obs['kmeans'],average_method='arithmetic')
#     ami_hc = adjusted_mutual_info_score(adata.obs['celltype'], adata.obs['hc'],average_method='arithmetic')
#     #homogeneity
#     homo_louvain = homogeneity_score(adata.obs['celltype'], adata.obs['louvain'])
#     homo_kmeans = homogeneity_score(adata.obs['celltype'], adata.obs['kmeans'])
#     homo_hc = homogeneity_score(adata.obs['celltype'], adata.obs['hc'])

#     df_metrics = pd.DataFrame(
#         columns=[
#             "ARI_Louvain",
#             "ARI_kmeans",
#             "ARI_HC",
#             "AMI_Louvain",
#             "AMI_kmeans",
#             "AMI_HC",
#             "Homogeneity_Louvain",
#             "Homogeneity_kmeans",
#             "Homogeneity_HC",
#         ]
#     )

#     df_metrics.loc[adata.uns["method"], ["ARI_Louvain", "ARI_kmeans", "ARI_HC"]] = [
#         ari_louvain,
#         ari_kmeans,
#         ari_hc,
#     ]
#     df_metrics.loc[adata.uns["method"], ["AMI_Louvain", "AMI_kmeans", "AMI_HC"]] = [
#         ami_louvain,
#         ami_kmeans,
#         ami_hc,
#     ]
#     df_metrics.loc[
#         adata.uns["method"],
#         ["Homogeneity_Louvain", "Homogeneity_kmeans", "Homogeneity_HC"],
#     ] = [homo_louvain, homo_kmeans, homo_hc]

#             print(
#                 "too large. Update max: {} -> {}".format(self._scan[1], self.current),
#                 end=". ",
#             )
#             print(
#                 "too small. Update min: {} -> {}".format(self._scan[0], self.current),
#                 end=". ",
#             )
#         print("New bounds: {}".format(self._scan), end=". ")
#         print("Current resolution: {}".format(self.current), end=". ")
# dcr = DynamicClusteringResolution(target=9, scan=(1, 11))
# dcr.current
# dcr.evaluate(12)

# ------------------------------------------------------------------------------------------------

# def update_resolution(scan, resolution, n_clusters_target, n_clusters_predicted):

#     if n_clusters_predicted:
#         if n_clusters_predicted > n_clusters_target


#     resolution = np.mean(scan)

#     return resolution, scan

# import time

# def message_tracer(abstract_func):

#     def message_wrapper(*args, **kwargs):

#         result = abstract_func(*args, **kwargs)

#         t_delta = tf - t0

#         if not result is None:
#             return result, t_delta
#         return t_delta

#     return message_wrapper

# class VerboseMessaging:
#     def __init__(self):
#         pass

#     @property
#     def too_large(self):

#     this_step = 0
#     this_min = float(range_min)
#     this_max = float(range_max)
#     while this_step < max_steps:
#         print('step ' + str(this_step))
# #         adata.X = adata.obsm['X_umap']
#         this_resolution = this_min + ((this_max-this_min)/2)
#         sc.tl.louvain(adata,resolution=this_resolution)
#         this_clusters = adata.obs['louvain'].nunique()
#         if this_clusters > n_cluster:
#             this_max = this_resolution
#         elif this_clusters < n_cluster:
#             this_min = this_resolution
#         else:
#             return(this_resolution, adata)
#         this_step += 1

# max_steps = 5
# scan = (1, 5)
# n_clusters = 9
# # -- loop: --------
# step = 0
# target = n_clusters
# while step < max_steps:
#     print("\nstep: {}".format(step))
#     current_min, current_max = scan[0], scan[1]
#     # start in the middle
#     res = current_min + ((current_max - current_min) / 2)
#     n_clusters = np.random.randint(low=5, high=15, size=1)
#     print(
#         "Target: {} | Created: {} | Resolution: {:.2f}".format(
#             target, n_clusters[0], res
#         )
#     )
#     if n_clusters != target:
#         # adjust res
#         if n_clusters > target:
#             print("updating min")
#             current_min = res
#         else:
#             print("updating max")
#             current_max = res
#         print("new min-max: {:.2f}-{:.2f}".format(current_min, current_max))

#     step += 1