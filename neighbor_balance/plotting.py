import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import EngFormatter
from matplotlib.gridspec import GridSpec
import numpy as np
import itertools
from .ice import get_distance_average, get_marginal
import logging
import matplotlib as mpl
import cooler
import pyBigWig
from scipy.ndimage import gaussian_filter1d, gaussian_filter
from scipy.signal import find_peaks

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


def load_bw(track, chrom, start, end, default=0):
    """
    Load a bigwig track and return the values in a given region.
    """
    bw = pyBigWig.open(track)
    if chrom not in bw.chroms():
        print(f'Chromosome {chrom} not found in {track}')
        return np.zeros(end-start) * np.nan
    vals = np.array(bw.values(chrom, start, end))
    vals[np.isnan(vals)] = default
    return vals


def bin_track(vals, step=200):
    """
    Coarse-grain a track into bins of size step.
    """
    return np.array([np.nanmean(vals[i:i+step]) for i in range(0, len(vals), step)])


def get_epigenetics(track, chrom, start, end, smoothing=0, bin=True):
    """
    Load a bigwig track and return the values in a given region.

    Parameters
    ----------
    track: str
        The path to the bigwig file.
    chrom: str
        The chromosome to load.
    start: int
        The start position of the region to load.
    end: int
        The end position of the region to load.
    smoothing: int
        The smoothing window size. If <= 0, no smoothing is applied.
    bin: bool
        Whether to bin the track into bins of size 200. This reduces the size of the array.

    Returns
    -------
    x: np.ndarray
        The x coordinates of the track.
    vals: np.ndarray
        The values of the track.
    """
    vals = load_bw(track, chrom, start, end)
    if smoothing > 0:
        vals = gaussian_filter1d(vals, smoothing)
    if bin:
        vals = bin_track(vals, 200)
        x = np.arange(start, end, 200)
    else:
        x = np.arange(start, end)
    return x, vals


def get_epigenetics_ylims(contact_maps, tracks):
    """
    Get the maximum values of the epigenetic tracks in the contact maps
    """
    ylims = {}
    for track_name, track in tracks.items():
        high = -float('inf')
        for contact_map in contact_maps:
            _, vals = get_epigenetics(track, contact_map.chrom, contact_map.start, contact_map.end)
            high = max(high, np.nanmax(vals))
        ylims[track_name] = (0, high)
    return ylims


class ContactMap:
    def __init__(self, contact_map, chrom, start, end, resolution):
        self.contact_map = contact_map
        self.chrom = chrom
        self.start = start
        self.end = end
        self.resolution = resolution

    @classmethod
    def from_npz(cls, npz_file):
        """
        Load a contact map from a npz file.

        This would typically be the output from computing an ICE balanced contact map for a specific genomic region.
        The npz file should contain two entries: 'array' for the contact map and 'metadata' for the region information.
        The metadata should include 'region' (str: chrom:start-end) and 'resolution' (int: 200).
        
        Parameters
        ----------
        npz_file: str
            The path to the npz file.

        Returns
        -------
        ContactMap
        """
        loaded = np.load(npz_file, allow_pickle=True)
        contact_map = loaded['array']
        metadata = loaded['metadata'].item()
        chrom, start, end = parse_region(metadata['region'])
        return ContactMap(contact_map, chrom, start, end, metadata['resolution'])

    @classmethod
    def from_cooler(cls, path: str, resolution: int, region: str, balance=True, min_alpha=-1,
                    cooler_internal_path='::/resolutions/') -> "ContactMap":
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
        ContactMap
        """
        print(f'{path}{cooler_internal_path}{resolution}')
        clr = cooler.Cooler(f'{path}{cooler_internal_path}{resolution}')
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

    def compare(self, other, self_name=None, other_name=None, zoom_start=None, zoom_end=None, vmin=1e-3, vmax=1, bw=0, density_max=70):
        if (self_name is None) != (other_name is None):
            raise ValueError('Both contact maps must have names or neither can have names.')
        if not self.are_comparable(other):
            raise ValueError('Contact maps are not comparable. They must have the same chromosome and start/end positions.')

        f, ax = plt.subplots(2, 9, figsize=(12, 3.75), sharex='col',
                            gridspec_kw={'hspace': 0.1, 'height_ratios': [10, 1], 'width_ratios': [20, 1, 3]*3, 'wspace': 0.05})

        def cleanup(i):
            if zoom_start is not None:
                ax[0, i].set_xlim(zoom_start, zoom_end)
                ax[0, i].set_ylim(zoom_end, zoom_start)
                ax[1, i].set_xlim(zoom_start, zoom_end)
            format_ticks(ax[1, i], y=False)
            ax[1, i].set_ylim(0, density_max)
            ax[1, i+1].axis('off')
            ax[0, i+2].axis('off')
            ax[1, i+2].axis('off')

        if bw > 0:
            self = self.copy()
            other = other.copy()
            self.contact_map = np.tril(gaussian_filter(self.contact_map, bw)) + np.triu(self.contact_map, k=1)
            other.contact_map = np.tril(gaussian_filter(other.contact_map, bw)) + np.triu(other.contact_map, k=1)
        
        i = 0
        im = other.plot_contact_map(ax=ax[0, i], vmin=vmin, vmax=vmax, colorbar=False)
        plt.colorbar(im, cax=ax[0, i+1], orientation='vertical', extend='min')
        ax[1, i].plot(other.x(), other.get_marginal(), c='gray')
        ax[1, i].set_ylabel('Contact density', rotation=0, ha='right', va='center')
        cleanup(i)
        if other_name is not None:
            ax[0, i].set_title(other_name)
    
        i = 3
        im = self.plot_contact_map(ax=ax[0, i], vmin=vmin, vmax=vmax, colorbar=False)
        plt.colorbar(im, cax=ax[0, i+1], orientation='vertical', extend='min')
        ax[1, i].plot(self.x(), self.get_marginal(), c='black')
        
        cleanup(i)
        ax[0, i].set_yticklabels([])
        ax[1, i].set_yticklabels([])
        if self_name is not None:
            ax[0, i].set_title(self_name)

        i = 6
        log_change = self.copy()
        log_change.contact_map = np.log2(self.contact_map / other.contact_map)
        im = log_change.plot_contact_map(cmap='coolwarm', vmin=-2, vmax=2, colorbar=False, ax=ax[0, i], log_norm=False)
        plt.colorbar(im, cax=ax[0, i+1], orientation='vertical')
        ax[1, i].plot(other.x(), other.get_marginal(), c='gray')
        ax[1, i].plot(self.x(), self.get_marginal(), c='black')
        cleanup(i)
        ax[0, i].set_yticklabels([])
        ax[1, i].set_yticklabels([])
        if self_name is not None and other_name is not None:
            ax[0, i].set_title(f'log2 {self_name} / {other_name}')

        return f, ax
    
    def epigenetics_plot(self, tracks, depth=None, vmin=1e-4, smoothing=200, bin=True, show_map=True, ylims=None,
                         contact_height=0.75, width=20, track_height=0.5, mean_density=None):
        if depth is None:
            depth = (self.end - self.start) / 3

        if show_map:
            height_ratios = [width * depth / (1.5 * (self.end - self.start))]
            height_ratios += [contact_height]
            height_ratios += [track_height]*len(tracks)
            f, axs = plt.subplots(len(height_ratios), 1, figsize=(width, sum(height_ratios)), sharex=True,
                                gridspec_kw={'height_ratios': height_ratios})
            self.plot_contact_map_horizontal(ax=axs[0], depth=depth, vmin=vmin, colorbar=False)
            contact_map_ax = axs[0]
            axs = axs[1:]
        else:
            height_ratios = [contact_height]
            height_ratios += [track_height]*len(tracks)
            f, axs = plt.subplots(len(height_ratios), 1, figsize=(width, sum(height_ratios)), sharex=True,
                                gridspec_kw={'height_ratios': height_ratios})
            contact_map_ax = None

        def format_ylabel(ax, name):
            # ax.set_ylabel(name, rotation=0, ha='left', labelpad=-10,
            #               bbox=dict(facecolor='white', edgecolor='none', alpha=0.5))
            ax.set_ylabel(name, rotation=0, ha='right', va='center')
            ax.yaxis.tick_right()
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            if ylims is not None:
                ax.set_ylim(*ylims[name])
                
        marginal = self.get_marginal()
        peaks, _ = find_peaks(gaussian_filter1d(-np.log2(marginal), 5), prominence=0.4)
        print(peaks)
        for peak in peaks:
            for ax in axs:
                ax.axvline(self.x()[peak], c='gray', lw=1, alpha=0.5)

        axs[0].set_xlim(self.start, self.end)
        format_ticks(axs[0], y=False)

        if mean_density is None:
            mean_density = np.nanmean(marginal)
        axs[0].axhline(mean_density, ls='--', c='grey')
        axs[0].set_yticks([mean_density])
        axs[0].plot(self.x(), marginal, c='black')
        format_ylabel(axs[0], 'Contact\ndensity')

        for i, (name, track) in enumerate(tracks.items()):
            ax = axs[i+1]
            x, vals = get_epigenetics(track, self.chrom, self.start, self.end, smoothing=smoothing, bin=bin)
            ax.fill_between(self.x(), np.zeros(vals.shape), vals, color='green')
            format_ylabel(ax, name)

            ax.set_ylim(0)
            ax.set_xlim((self.start, self.end))
            ax.get_yticklabels()[0].set_visible(False)
        return f, axs, contact_map_ax


    def plot_contact_map_and_marginal(self, vmax=1, vmin=0.00001):
        margin = 0.5
        sep = 0.2
        base = 0.4
        h = 22 * base + 2 * margin + sep
        w = 21 * base + 3 * margin + sep
        left_right = margin / w
        top_bottom = margin / h
        f = plt.figure(figsize=(w, h), dpi=200)
        gs = GridSpec(2, 2, height_ratios=[20, 2], width_ratios=[20, 1],
                    left=2*left_right, right=1-left_right, wspace=sep/w,
                    top=1-top_bottom, bottom=top_bottom, hspace=sep/h)

        image_ax = f.add_subplot(gs[0, 0])
        cbar_ax = f.add_subplot(gs[0, 1].subgridspec(3, 1, height_ratios=[0.25, 1, 0.25])[1])
        plot_ax = f.add_subplot(gs[1, 0], sharex=image_ax)

        im = self.plot_contact_map(ax=image_ax, vmax=vmax, vmin=vmin, colorbar=False)
        plt.colorbar(im, cax=cbar_ax, extend='min')
        image_ax.tick_params(labelbottom=False)

        plot_ax.plot(self.x(), self.get_marginal(), c='black')
        plot_ax.tick_params(axis='x', rotation=20)
