import matplotlib.pyplot as plt
import pandas as pd
from importlib.resources import files

# %% ----
# CLASSES
# -------


class Ideogram:
    """
    Ideogram object
    """

    def __init__(
        self,
        chroms,
        n_tracks_above=0,
        n_tracks_below=0,
        chrom_height=0.6,
        track_height=0.2,
        track_gap=0.05,
        chrom_track_gap=0.1,
        row_gap=0.3,
    ):
        self.chroms = chroms
        self.n_tracks_above = n_tracks_above
        self.n_tracks_below = n_tracks_below
        self.chrom_height = chrom_height
        self.track_height = track_height
        self.track_gap = track_gap
        self.chrom_track_gap = chrom_track_gap
        self.row_gap = row_gap

        self.above_span = self._span(self.n_tracks_above)
        self.below_span = self._span(self.n_tracks_below)

        self.row_height = self.above_span + self.chrom_height + self.below_span
        self.row_stride = self.row_height + self.row_gap

        self.chrom_y_center = {
            chrom: i * self.row_stride + self.above_span + self.chrom_height / 2
            for i, chrom in enumerate(self.chroms)
        }

    def _span(self, n_tracks):
        """
        Shared formula for computing the vertical space occupied by a stack
        of tracks (used for both the 'above' and 'below' stacks).
        """
        return (
            n_tracks * self.track_height
            + max(n_tracks - 1, 0) * self.track_gap
            + (self.chrom_track_gap if n_tracks else 0)
        )

    def chrom_yrange(self, chrom):
        """
        Function that gets the yranges parameter for broken_barh (ypos + yextent)
        """

        yc = self.chrom_y_center[chrom]
        return (yc - self.chrom_height / 2, self.chrom_height)

    def track_yrange(self, chrom, position, index=0):
        """
        Function that gets the yranges param for a given track
        """

        yc = self.chrom_y_center[chrom]
        if position == "above":
            y_top = yc - self.chrom_height / 2 - self.chrom_track_gap
            y_top -= index * (self.track_height + self.track_gap)
            return (y_top - self.track_height, self.track_height)
        elif position == "below":
            y_bottom = yc + self.chrom_height / 2 + self.chrom_track_gap
            y_bottom += index * (self.track_height + self.track_gap)
            return (y_bottom, self.track_height)
        else:
            raise ValueError("position of track must be 'above' or 'below'")

    @property
    def total_height(self):
        return len(self.chroms) * self.row_stride - self.row_gap

    @property
    def yticks_and_labels(self):
        return list(self.chrom_y_center.values()), list(self.chrom_y_center.keys())


# %% ------
# FUNCTIONS
# ---------


def load_cytobands(genome="hg38"):
    data_path = files("karpyoplot.data") / f"{genome}_cytobands.bed"
    return pd.read_csv(
        data_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "stop", "band", "stain"],
    )


def _valid_chrom_groups(df, layout, chrom_col="chr"):
    """
    Helper function that iterates over (chrom, subdf) groups, skipping any
    chromosome not present in the layout.
    """

    for chrom, subdf in df.groupby(chrom_col):
        if chrom in layout.chrom_y_center:
            yield chrom, subdf


def plot_ideogram(ax, cyto, layout, chrom_col="chr"):
    """
    Function that prepares the fig
    """

    for chrom, subdf in _valid_chrom_groups(cyto, layout, chrom_col=chrom_col):
        subdf = subdf.sort_values(by="start")
        xranges = list(zip(subdf["start"], subdf["stop"] - subdf["start"]))
        colors = subdf["color"].tolist()
        ax.broken_barh(
            xranges,
            layout.chrom_yrange(chrom),
            facecolors=colors,
            edgecolors="black",
            linewidths=0.2,
        )


def plot_region(
    ax,
    layout,
    df,
    position,
    index=0,
    color="blue",
    linewidth=0.2,
    label=None,
    fontsize=6,
    rasterized=False,
    chrom_col="chr",
):
    """
    Function that can add regions as a track. Needs bed like format with chr, start, stop.
    Can modify the "chr" col with `chrom_col`
    """

    for chrom, subdf in _valid_chrom_groups(df, layout, chrom_col=chrom_col):
        yrange = layout.track_yrange(chrom, position=position, index=index)
        xranges = list(zip(subdf["start"], subdf["stop"] - subdf["start"]))
        ax.broken_barh(
            xranges,
            yrange,
            facecolors=color,
            linewidths=linewidth,
            edgecolor=color,
            rasterized=rasterized,
        )

    if label:
        for chrom in layout.chroms:
            y_bottom, height = layout.track_yrange(
                chrom, position=position, index=index
            )
            ax.annotate(
                label,
                xy=(0, y_bottom + height / 2),
                xytext=(-4, 0),
                textcoords="offset points",
                ha="right",
                va="center",
                fontsize=fontsize,
            )


def _rgb_to_hex(r, g, b):
    """
    Helper function that converts rgb code to hex code
    """

    return "#%02x%02x%02x" % (r, g, b)


def plot_ucsc_track(
    ax,
    layout,
    df,
    position,
    index=0,
    color="black",
    color_col=None,
    linewidth=0.2,
    label=None,
    fontsize=6,
    rasterized=False,
    chrom_col="chr",
):
    """
    Function very similar to plot_region(), but accepts bed file directly from bigBedToBed
    """

    df = df.copy()

    if color_col:
        temp_rgb_df = df[color_col].str.split(",", expand=True).astype(int)
        df["hex_color"] = temp_rgb_df.apply(lambda r: _rgb_to_hex(*r), axis=1)

    # print(df.head())

    for chrom, subdf in _valid_chrom_groups(df, layout, chrom_col=chrom_col):
        yrange = layout.track_yrange(chrom, position=position, index=index)
        xranges = list(zip(subdf["start"], subdf["stop"] - subdf["start"]))

        if color_col:
            colors = subdf["hex_color"].tolist()
        else:
            colors = color

        ax.broken_barh(
            xranges,
            yrange,
            facecolors=colors,
            linewidths=linewidth,
            edgecolor=colors,
            rasterized=rasterized,
        )

    if label:
        for chrom in layout.chroms:
            y_bottom, height = layout.track_yrange(
                chrom, position=position, index=index
            )
            ax.annotate(
                label,
                xy=(0, y_bottom + height / 2),
                xytext=(-4, 0),
                textcoords="offset points",
                ha="right",
                va="center",
                fontsize=fontsize,
            )


def plot_marker(
    ax,
    layout,
    position,
    index,
    color="black",
    linewidth=1,
    fontsize=6,
    df=None,
    name_col="gene_name",
    strand_col="strand",
    chrom_col="chr",
    label=None,
    chrom=None,
    pos=None,
    strand=None,
    rotation=0,
    arrow_length=2_000_000,
    arrow_height=None,
    gap_override=None,
    fontstyle=None,
):
    """
    Function to plot markers as small arrows + labels (genes, etc.) on a track.

    Two modes:
      - df given: plots one marker+label per row (needs chrom_col, start, optionally stop,
      strand_col, name_col). if stop present, marker is placed between start and stop
      - df=None: plots a single marker+label at (chrom, pos, strand), labeled with `label`

    gap_override: if given, overrides layout.chrom_track_gap for this marker's base point,
    letting the arrow start closer to (or on, with 0) the chromosome bar rather than at
    the track slot's normal near edge.

    arrow_height: if given, sets the length of the arrow's vertical segment directly
    (in data/y-axis units), independent of track_height.

    arrow_length: controls the length of the arrow's pointy end.

    arrow_height and arrow_length are useful when it comes to arranging gene labels and other
    aesthetics
    """

    if df is None and (chrom is None or pos is None or strand is None):
        raise ValueError("Must provide either `df`, or `chrom`, `pos`, and `strand`.")

    def _draw_one(c, x, s, text):
        if c not in layout.chrom_y_center:
            return
        y_bottom, height = layout.track_yrange(c, position=position, index=index)
        y_base = y_bottom + height if position == "above" else y_bottom
        y_tip = y_bottom if position == "above" else y_bottom + height

        if gap_override is not None:
            chrom_bottom, chrom_height = layout.chrom_yrange(c)
            if position == "above":
                y_base = chrom_bottom - gap_override
            else:
                y_base = chrom_bottom + chrom_height + gap_override

        if arrow_height is not None:
            y_tip = (
                y_base - arrow_height if position == "above" else y_base + arrow_height
            )

        angleA = -90 if position == "above" else 90
        dx = arrow_length if s == "+" else -arrow_length
        ha = "left" if s == "+" else "right"

        ax.annotate(
            "",
            xy=(x + dx, y_tip),
            xytext=(x, y_base),
            arrowprops=dict(
                arrowstyle="-|>",
                color=color,
                linewidth=linewidth,
                connectionstyle=f"angle,angleA={angleA},angleB=0,rad=0",
                shrinkA=0,
                shrinkB=0,
            ),
        )

        ax.annotate(
            text,
            xy=(x + dx, y_tip),
            xytext=(2 if s == "+" else -2, 0),
            textcoords="offset points",
            ha=ha,
            va="center",
            fontsize=fontsize,
            rotation=rotation,
            fontstyle=fontstyle,
        )

    if df is not None:
        for c, subdf in df.groupby(chrom_col):
            positions = (
                (subdf["start"] + subdf["stop"]) / 2
                if "stop" in subdf.columns
                else subdf["start"]
            )
            for x, s, name in zip(positions, subdf[strand_col], subdf[name_col]):
                _draw_one(c, x, s, name)
    else:
        _draw_one(chrom, pos, strand, label)


# %%
