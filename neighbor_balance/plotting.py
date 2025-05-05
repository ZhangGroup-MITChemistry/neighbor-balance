import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import EngFormatter
import numpy as np
import itertools
from .ice import get_distance_average, get_marginal
import logging
import matplotlib as mpl
import cooler

# cooltools.lib.plotting is currently not compatible with matplotlib > 3.8. All we need from it is the fall color map,
# so I've copied the definition here. cooltools will probably be fixed in the next release, so at that point this
# can be removed and the import uncommented.

# import cooltools.lib.plotting  # need to import this to get the fall color map.
pal = np.array(
    (
        (255, 255, 255),
        (255, 255, 204),
        (255, 237, 160),
        (254, 217, 118),
        (254, 178, 76),
        (253, 141, 60),
        (252, 78, 42),
        (227, 26, 28),
        (189, 0, 38),
        (128, 0, 38),
        (0, 0, 0),
    )
)
mpl.colormaps.register(mpl.colors.LinearSegmentedColormap.from_list('fall', pal / 255.0, 256))

def apply_matplotlib_style():
    logging.basicConfig(level=logging.CRITICAL)
    logging.getLogger("fontTools.subset").setLevel(logging.CRITICAL)
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['ytick.labelsize'] = 8
    mpl.rcParams['xtick.labelsize'] = 8
    mpl.rcParams['legend.fontsize'] = 8
    mpl.rcParams['figure.dpi'] = 150
    mpl.rcParams['savefig.transparent'] = True
    mpl.rcParams['savefig.bbox'] = 'tight'


def parse_region(region: str) -> tuple[str, int, int]:
    chrom, rest = region.split(':')
    start, end = rest.split('-')
    start = int(start.replace(',', ''))
    end = int(end.replace(',', ''))
    return chrom, start, end


def format_ticks(ax, x=True, y=True, rotate=True):
    bp_formatter = EngFormatter('b')

    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x', rotation=20)


def _pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
    # From cooltools documentation.
    start_pos_vector = [start+resolution*i for i in range(len(matrix_c)+1)]
    n = matrix_c.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    matrix_a = np.dot(np.array([(i[1], i[0])
                                for i in itertools.product(start_pos_vector[::-1],
                                                           start_pos_vector)]), t)
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)
    im = ax.pcolormesh(x, y, np.flipud(matrix_c), *args, **kwargs)
    im.set_rasterized(True)
    return im


def label_comparison(ax, name1, name2, fontsize=16, **text_kwargs):
    ax.text(0.99, 0.99, name1, transform=ax.transAxes, ha='right', va='top', fontsize=fontsize, **text_kwargs)
    ax.text(0.01, 0.01, name2, transform=ax.transAxes, ha='left', va='bottom', fontsize=fontsize, **text_kwargs)


class ContactMap:
    def __init__(self, contact_map, chrom, start, end, resolution):
        self.contact_map = contact_map
        self.chrom = chrom
        self.start = start
        self.end = end
        self.resolution = resolution

    @classmethod
    def from_cooler(cls, path: str, resolution: int, region: str, balance=True, min_alpha=-1) -> "ContactMap":
        """
        Load a contact map from a cooler file.

        Parameters
        ----------
        path: str
            The path to the cooler file.
        resolution: int
            The resolution of the contact map.
        region: str
            The region of the contact map.
        balance: bool
            Whether to load the ice balanced contact map or raw counts.
        min_alpha: int
            Minimum value balancing weight, used as a proxy for row coverage.
            This option is deprecated and should not be used. If there are very
            high variance rows with low coverage, go back and run ICE balancing with
            a smaller `mad_max` value or filter regions not covered with capture probes.

        Returns
        -------
        np.ndarray
            The contact map.
        """
        clr = cooler.Cooler(f'{path}::resolutions/' + str(resolution))
        contact_map = clr.matrix(balance=balance).fetch(region)
        if min_alpha > 0:
            logging.warning('Use of min_alpha is deprecated. Filter rows the right way.')
            i, j = clr.extent(region)
            h5 = clr.open()
            alpha = 1 / h5['bins']['weight'][i:j]
            mask = alpha < min_alpha
            contact_map[mask, :] = np.nan
            contact_map[:, mask] = np.nan

        chrom, start, end = parse_region(region)
        return ContactMap(contact_map, chrom, start, end, resolution)

    def x(self):
        return np.arange(self.start, self.end, self.resolution) + self.resolution / 2

    def plot_contact_map(self, ax=None, colorbar=True, vmax=1, vmin=1e-5,
                         log_norm=True, cmap='fall', extend='min'):

        if ax is None:
            f, ax = plt.subplots()

        if log_norm:
            norm = LogNorm(vmax=vmax, vmin=vmin)
        else:
            norm = Normalize(vmax=vmax, vmin=vmin)

        im = ax.matshow(
            self.contact_map,
            cmap=cmap,
            norm=norm,
            extent=(self.start, self.end, self.end, self.start),
            rasterized=True,
        )

        if colorbar:
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, extend=extend)
        format_ticks(ax)
        return im

    def plot_contact_map_horizontal(self, ax=None, colorbar=True, vmax=1, vmin=1e-5,
                                    log_norm=True, cmap='fall', extend='min', depth=None):

        if depth is None:
            depth = self.end - self.start

        if ax is None:
            f, ax = plt.subplots()

        if log_norm:
            norm = LogNorm(vmax=vmax, vmin=vmin)
        else:
            norm = Normalize(vmax=vmax, vmin=vmin)

        im = _pcolormesh_45deg(ax, self.contact_map, start=self.start, resolution=self.resolution, cmap=cmap, norm=norm)

        if colorbar:
            plt.colorbar(im, ax=ax, extend=extend)
        format_ticks(ax)
        ax.set_aspect(0.5)
        ax.set_ylim(0, depth)
        return im

    def plot_distance_average(self, ax=None, show_self=False, **plot_kwargs):
        if ax is None:
            f, ax = plt.subplots()
        y = get_distance_average(self.contact_map, show_self=show_self)
        shift = 0 if show_self else 1
        ax.plot(self.resolution * np.arange(shift, len(y) + shift), y, **plot_kwargs)
        ax.set_xscale("log")
        ax.set_yscale("log")
        format_ticks(ax, y=False)

    def are_comparable(self, other):
        return (self.chrom == other.chrom) and (self.start == other.start) and (self.end == other.end)

    def get_merged_map(self, other, use_main_diagonal=True):
        assert self.are_comparable(other)
        merged_map = np.triu(self.contact_map, k=0 if use_main_diagonal else 1) + np.triu(other.contact_map, k=1).T
        return ContactMap(merged_map, self.chrom, self.start, self.end, self.resolution)

    def copy(self):
        return ContactMap(self.contact_map.copy(), self.chrom, self.start, self.end, self.resolution)

    def select_region(self, start, end):
        assert (start - self.start) % self.resolution == 0
        assert (end - self.end) % self.resolution == 0
        assert start >= self.start and end >= self.end

        start_i = (start - self.start) // self.resolution
        end_i = (end - self.start) // self.resolution
        return ContactMap(self.contact_map[start_i:end_i, start_i:end_i], self.chrom, start, end, self.resolution)

    def get_marginal(self, k=1, correct_for_flanks=False):
        return get_marginal(self.contact_map, k=k, correct_for_flanks=correct_for_flanks)


# def plot_contact_map_and_marginal(contact_map: np.ndarray, region: str, vmax=1, vmin=0.00001):
#     margin = 0.5
#     sep = 0.2
#     base = 0.4
#     h = 22 * base + 2 * margin + sep
#     w = 21 * base + 3 * margin + sep
#     left_right = margin / w
#     top_bottom = margin / h
#     f = plt.figure(figsize=(w, h), dpi=200)
#     gs = GridSpec(2, 2, height_ratios=[20, 2], width_ratios=[20, 1],
#                   left=2*left_right, right=1-left_right, wspace=sep/w,
#                   top=1-top_bottom, bottom=top_bottom, hspace=sep/h)
#
#     image_ax = f.add_subplot(gs[0, 0])
#     cbar_ax = f.add_subplot(gs[0, 1].subgridspec(3, 1, height_ratios=[0.25, 1, 0.25])[1])
#     plot_ax = f.add_subplot(gs[1, 0], sharex=image_ax)
#
#     im = plot_contact_map(contact_map, region, ax=image_ax, vmax=vmax, vmin=vmin, colorbar=False)
#     plt.colorbar(im, cax=cbar_ax, extend='min')
#     image_ax.tick_params(labelbottom=False)
#
#     chrom, start, end = parse_region(region)
#     resolution = (end - start) // contact_map.shape[0]
#     plot_ax.plot(range(start, end, resolution), get_marginal(contact_map))
#     plot_ax.tick_params(axis='x', rotation=20)
