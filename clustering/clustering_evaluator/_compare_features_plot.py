
from ._custom_plot_layout import CustomPlotLayout

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import vinplots

xlims = {
    "10xpbmc_5k": (0.3, 0.7),
    "BMsim_clean": (0.9, 1.1),
    "BMsim_noisyp4": (0.9, 1.1),
    "buenrostro2018": (0.3, 0.8),
    "sciatac_subset": (0.2, 0.7),
}


def group_subset(group_df, subset):
    return (
        pd.concat([group_df["features"], group_df.filter(regex=subset)], axis=1)
        .set_index("features")
        .T
    )


def header_formatting(
    ax, x=0.5, y=0.4, s=None, fc="#ECECEC", fontsize=8, ha="center", rotation=0
):
    ax.set_facecolor(fc)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 1.1)
    ax.text(x=x, y=y, s=s, fontsize=fontsize, ha=ha, rotation=rotation)


def _vertical_headers(nplots, ncols):
    return (np.arange(0, nplots + ncols)[::ncols] - 1)[2:].tolist()


def connector(ax, feature_1, feature_2, zorder=4, c="k", lw=1, **kwargs):
    for i in range(3):
        x1, x2 = (feature_1[i], feature_2[i]), (i, i)
        ax.plot(x1, x2, zorder=zorder, c=c, lw=lw, **kwargs)
        
def extract_labels(df, split_on="_", option_a=True):
    if option_a:
        return [split_on.join(val.split(split_on)[1:]) for val in df.index.tolist()]
    return [val.split(split_on)[0] for val in df.index.tolist()]


def plot_feature_comparison(
    ax,
    group_df_subset,
    feature_keys=["peaks", "peaks_sequences"],
    colors=["crimson", "dodgerblue"],
    split_on="_",
    s=60,
    ec="None",
):

    labels = extract_labels(group_df_subset, split_on=split_on)

    for n, key in enumerate(feature_keys):
        x = group_df_subset[key]
        y = range(len(group_df_subset[key]))
        c = ["w", colors[n]]
        z = [5, 6]
        for i in range(2):
            ax.scatter(x, y, s=s, ec=ec, c=c[i], alpha=1, zorder=z[i])

    connector(ax, group_df_subset[feature_keys[0]], group_df_subset[feature_keys[1]])
    
def format_data_ax(
    ax,
    dataset,
    clustering_methods,
    grid_color="#ECECEC",
    use_xticks=False,
    use_yticks=False,
):

    ax.grid(True, color=grid_color, zorder=0)

    # -- x axis: ------------------------------------------------------------------------
    xl1, xl2 = xlims[dataset][0], xlims[dataset][1]
    edge_space = (xl2 - xl1) / 10

    xticks = np.round(np.linspace(xl1, xl2, 5), 2)
    if use_xticks:
        ax.set_xticks(xticks, ["{:.2f}".format(xt) for xt in xticks])
    else:
        ax.set_xticks(xticks, [""] * len(xticks))
        ax.tick_params(axis="x", which="major", length=0.01)

    ax.set_xlim(xl1 - edge_space, xl2 + edge_space)
    ax.tick_params(axis="x", which="major", labelsize=6)

    # -- y axis: ------------------------------------------------------------------------
    ax.set_ylim(-0.5, len(clustering_methods) - 0.5)
    if use_yticks:
        ax.set_yticks(range(len(clustering_methods)), clustering_methods)
        ax.tick_params(axis="y", which="major", labelsize=8)
    else:
        ax.set_yticks(range(len(clustering_methods)), [""] * len(clustering_methods))
        ax.tick_params(axis="y", which="major", length=0.01)
        
class CompareFeaturePlot:
    def __init__(
        self,
        score_df,
        dataset_key="dataset",
        figsize=1.2,
        wspace=0.05,
        hspace=0.05,
        header_size_ratio=0.25,
        metrics=["ARI", "AMI", "Homogeneity"],
        clustering_methods=["louvain", "kmeans", "h_clust"],
    ):
        cp = CustomPlotLayout(
            figsize=np.array([6, 3]) * figsize,
            upperbound=0.90,
            rightbound=0.95,
            ncols=5,
            nrows=3,
            hspace=hspace,
            wspace=wspace,
        )
        self.fig, self.axes = cp()

        self.__parse__(locals())
        self.grouped = self.score_df.groupby(dataset_key)

        self.datasets = list(self.grouped.groups.keys())
        self.n_datasets = self.grouped.ngroups
        self.n_metrics = len(self.metrics)
        self.ncols = self.n_datasets + 1
        self.nrows = self.n_metrics + 1
        self.nplots = self.ncols * self.nrows
        width_ratios = [1] * self.n_datasets + [header_size_ratio]
        height_ratios = [header_size_ratio] + [1] * self.n_metrics

    def __parse__(self, kwargs, ignore=["self"]):

        for key, val in kwargs.items():
            if not key in ignore:
                setattr(self, key, val)

    def format_col_headers(self):

        for i, ax in enumerate(self.axes["Columns"]):
            ax.set_facecolor("#ECECEC")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlim(-0.1, 1.1)
            ax.set_ylim(-0.1, 1.1)
            ax.text(x=0.5, y=0.25, s=self.datasets[i], fontsize=8, ha="center")

    def format_row_headers(self):

        for i, ax in enumerate(self.axes["Rows"]):
            ax.set_facecolor("#ECECEC")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlim(-0.1, 1.1)
            ax.set_ylim(-0.1, 1.1)
            if i == 2:
                y_ = 0.08
            else:
                y_ = 0.4

            ax.text(
                x=0.5,
                y=y_,
                s=self.metrics[i],
                fontsize=8,
                ha="center",
                rotation=90,
            )

    def forward(self, ax, group_df, dataset, metric, xticks, yticks):
        plot_feature_comparison(ax=ax, group_df_subset=group_subset(group_df, metric))
        format_data_ax(
            ax, dataset, self.clustering_methods, use_xticks=xticks, use_yticks=yticks
        )

    def __call__(self):

        self.format_col_headers()
        self.format_row_headers()
        for i, metric in enumerate(self.metrics):
            for j, (dataset, group_df) in enumerate(self.grouped):
                self.forward(
                    self.axes["Data"][(i, j)],
                    group_df,
                    dataset,
                    metric,
                    xticks=(i == 2),
                    yticks=(j == 0),
                )
                
                
def compare_features_plot(
    score_df, save=False, savename="clustering_metrics.feature_comparison.svg", **kwargs
):

    CFP = CompareFeaturePlot(score_df, **kwargs)
    CFP()

    if save:
        plt.savefig(savename)