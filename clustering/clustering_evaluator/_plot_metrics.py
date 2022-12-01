
import vinplots
import numpy as np
import matplotlib.pyplot as plt


class MultiMarkerScatter:
    MarkerDict = {"h_clust": "o", "kmeans": "^", "louvain": "s"}
    _colors = {
        "h_clust": "blue",
        "kmeans": "green",
        "louvain": "orange",
    }

    def __call__(self, ax, df, y, include_label=False):
        for idx in df.index:
            x = df.loc[idx]
            m = self.MarkerDict[idx]
            c = self._colors[idx]
            if include_label:
                l = idx
            else:
                l = None
            img = ax.scatter(x, y, marker=m, label=l, c=c, zorder=2)

        return img, ax


class MetricsFigure:
    metrics = ["ARI", "AMI", "Homogeneity"]

    MMScatter = MultiMarkerScatter()
    _counter = 0

    def __init__(self):
        self.fig = vinplots.Plot()
        self.fig.construct(
            nplots=6,
            ncols=3,
            height_ratios=[0.08, 1.0],
            figsize_height=0.6,
            figsize_width=0.5,
            wspace=0.3,
            hspace=0,
        )
        self.axes = self.fig.linearize()
        self.__format__()

    def __format__(self):

        # -- format title-box headers: ---------------------------------------------------
        for i_1 in range(3):
            self.axes[i_1].set_facecolor("#ECECEC")
            self.axes[i_1].set_xticks([])
            self.axes[i_1].set_yticks([])
            self.axes[i_1].set_xlim(-0.1, 1.1)
            self.axes[i_1].set_ylim(-0.1, 1.1)
            self.axes[i_1].text(
                x=0.5, y=0.25, s=self.metrics[i_1], fontsize=10, ha="center"
            )

        # -- format data plots: ----------------------------------------------------------
        for i_2 in range(3, 6):
            self.axes[i_2].set_xlim(-0.1, 1.1)
            self.axes[i_2].grid(True, color="#ECECEC", zorder=0)
            self.axes[i_2].tick_params(axis="both", which="major", labelsize=6)
            self.axes[i_2].set_xticks(np.linspace(0, 1, 5))

    def plot(self, adata, key="clustering_scores", add_label=False):

        df = adata.uns[key]
        include_label = False
        for n, col in enumerate(df.columns):
            if (n == 2) and (add_label):
                include_label = True
            img, ax = self.MMScatter(
                self.axes[n + 3], df[col], y=self._counter, include_label=include_label
            )

        if include_label:
            ax.legend(edgecolor="w", loc=(1.1, 0.9))
        self._counter += 1

    # reformat y-axis limits after plotting
    # if wrapped in a func, can do this on the fly after all methods have been plotted
    # # plot.axes[3].set_ylim(-1, 2)
    
    
def plot_metrics(DataDict, add_label = False, save=False, savename="clustering_metrics.feature_comparison.svg"):
    
    plot = MetricsFigure()
    
    data_keys = np.sort(list(DataDict.keys())).tolist()
    n_datasets = len(data_keys)

    for i, name in enumerate(data_keys):
        adata = DataDict[name]
        if i == (n_datasets - 1):
            add_label = True
        plot.plot(adata, add_label=add_label)
    xt = plot.axes[3].set_yticks(range(n_datasets), data_keys)
    
    if save:
        plt.savefig(savename)