import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import EngFormatter
from matplotlib.gridspec import GridSpec
import numpy as np
import itertools
from .ice import get_marginal, get_distance_average
import logging
import matplotlib as mpl
# import cooltools.lib.plotting


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


def label_comparison(ax, name1, name2):
    ax.text(0.99, 0.99, name1, transform=ax.transAxes, ha='right', va='top', fontsize=16)
    ax.text(0.01, 0.01, name2, transform=ax.transAxes, ha='left', va='bottom', fontsize=16)


def compare_on_diagonals(map1, map2, name1, name2, region, ax=None,
                         vmin=-0.3, vmax=0.3, cmap='coolwarm', main_diag_from_map1=True):
    if ax is None:
        f, ax = plt.subplots()
    merged = np.triu(map1, k=0 if main_diag_from_map1 else 1) + np.triu(map2, k=1).T
    plot_pairwise(merged, region, cmap=cmap, vmin=vmin, vmax=vmax, ax=ax)
    label_comparison(ax, name1, name2)


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


def plot_pairwise_horizontal(contact_map: np.ndarray, region: str, ax=None, colorbar=True,
                             vmax=None, vmin=None, log_norm=True, cmap='fall', extend='neither', depth=None):
    chrom, start, end = parse_region(region)
    resolution = (end - start) // contact_map.shape[0]
    if depth is None:
        depth = end - start

    if ax is None:
        f, ax = plt.subplots()

    if log_norm:
        norm = LogNorm(vmax=vmax, vmin=vmin)
    else:
        norm = Normalize(vmax=vmax, vmin=vmin)

    im = _pcolormesh_45deg(ax, contact_map, start=start, resolution=resolution, cmap=cmap, norm=norm)

    if colorbar:
        plt.colorbar(im, ax=ax, extend=extend)
    format_ticks(ax)
    ax.set_aspect(0.5)
    ax.set_ylim(0, depth)
    return im


def plot_pairwise(x: np.ndarray, region: str, ax=None, colorbar=True, vmax=None, vmin=None,
                  log_norm=False, cmap='hot_r', extend='neither'):
    chrom, start, end = parse_region(region)

    if ax is None:
        f, ax = plt.subplots()

    if log_norm:
        norm = LogNorm(vmax=vmax, vmin=vmin)
    else:
        norm = Normalize(vmax=vmax, vmin=vmin)

    im = ax.matshow(
        x,
        cmap=cmap,
        norm=norm,
        extent=(start, end, end, start),
        rasterized=True,
    )

    if colorbar:
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, extend=extend)
    format_ticks(ax)
    return im


def plot_contact_map(contact_map: np.ndarray, region: str, ax=None, colorbar=True, vmax=1, vmin=0.00001, cmap='fall'):
    return plot_pairwise(contact_map, region, ax=ax, colorbar=colorbar, vmax=vmax, vmin=vmin,
                         log_norm=True, cmap=cmap, extend='min')


def plot_contact_map_and_marginal(contact_map: np.ndarray, region: str, vmax=1, vmin=0.00001):
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

    im = plot_contact_map(contact_map, region, ax=image_ax, vmax=vmax, vmin=vmin, colorbar=False)
    plt.colorbar(im, cax=cbar_ax, extend='min')
    image_ax.tick_params(labelbottom=False)

    chrom, start, end = parse_region(region)
    resolution = (end - start) // contact_map.shape[0]
    plot_ax.plot(range(start, end, resolution), get_marginal(contact_map))
    plot_ax.tick_params(axis='x', rotation=20)


def plot_distance_average(contact_map, label, ax=None, show_self=False, color=None, bp_per_bin=200):
    if ax is None:
        f, ax = plt.subplots()
    y = get_distance_average(contact_map, show_self=show_self)
    shift = 0 if show_self else 1
    ax.plot(bp_per_bin * np.arange(shift, len(y)+shift), y, c=color, label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    format_ticks(ax, y=False)
