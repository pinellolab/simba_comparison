

# -- import packages: ------------
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import vinplots


# -- module class: ---------------
class CustomPlotLayout:
    def __parse__(self, kwargs, ignore=["self", "fig"]):
        for k, v in kwargs.items():
            if not k in ignore:
                setattr(self, "_{}".format(k), v)

    def __init__(
        self,
        figsize=(8, 3),
        upperbound=0.90,
        rightbound=0.95,
        ncols=5,
        nrows=3,
        hspace=0.3,
        wspace=0.3,
    ):

        self.fig = plt.figure(figsize=figsize)
        self.__parse__(locals())

    @property
    def header_row(self):

        self._header_row = plt.GridSpec(
            nrows=1,
            ncols=self._ncols,
            bottom=self._upperbound,
            top=1,
            left=0,
            right=self._rightbound,
            wspace=self._wspace,
            width_ratios=[1] * self._ncols,
        )

        return [
            self.fig.add_subplot(self._header_row[0, i])
            for i in range(self._header_row.ncols)
        ]

    @property
    def header_col(self):

        self._header_col = plt.GridSpec(
            nrows=self._nrows,
            ncols=1,
            bottom=0,
            top=self._upperbound,
            left=self._rightbound,
            right=1,
            hspace=self._hspace,
            height_ratios=[1] * self._nrows,
        )

        return [
            self.fig.add_subplot(self._header_col[i, 0])
            for i in range(self._header_col.nrows)
        ]

    @property
    def results(self):
        self._results = plt.GridSpec(
            nrows=self._nrows,
            ncols=self._ncols,
            bottom=0,
            top=self._upperbound,
            right=self._rightbound,
            left=0,
            hspace=self._hspace,
            wspace=self._wspace,
        )

        DataPlots = {}
        for i in range(self._results.nrows):
            for j in range(self._results.ncols):
                DataPlots[(i, j)] = self.fig.add_subplot(self._results[i, j])

        return DataPlots

    def __call__(self):
        self.Axes = {}
        self.Axes["Columns"] = self.header_row
        self.Axes["Rows"] = self.header_col
        self.Axes["Data"] = self.results

        return self.fig, self.Axes