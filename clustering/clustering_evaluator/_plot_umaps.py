
import vinplots


class UMAPFigure:
    def __init__(self, adata, nplots=12, ncols=4, **kwargs):

        self.adata = adata
        self.cols = adata.obs.columns.tolist()
        nplots = len(self.cols)
        self.scores = self.adata.uns["clustering_scores"]
        self.fig, self.axes = vinplots.quick_plot(
            nplots=nplots,
            ncols=ncols,
            figsize=0.8,
            wspace=0.2,
            hspace=0.2,
            spines_to_delete=["top", "right", "bottom", "left"],
            rm_ticks=True,
        )

    def _ax_formatting(self, ax, title=None, title_y=1.02):

        ax.set_title(title, y=title_y, fontsize=10)
        ax.set_xlabel("UMAP-1", fontsize=8)
        ax.set_ylabel("UMAP-2", fontsize=8)

    def plot_ax(self, ax, title=None, title_y=1.02):

        self._ax_formatting(ax, title=title, title_y=title_y)
        grouped = self.adata.obs.groupby(title)
        for group, group_df in grouped:
            group_adata = self.adata[group_df.index]
            xu = group_adata.obsm["X_umap"]
            ax.scatter(xu[:, 0], xu[:, 1], label=group, s=2)

    def __call__(self):

        for n, ax in enumerate(self.axes):
            self.plot_ax(ax, title=self.cols[n], title_y=1.02)