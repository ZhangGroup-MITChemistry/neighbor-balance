import click
import gzip
import os
import glob
import sys
from typing import Dict, Tuple, Optional, List

import numpy as np
import matplotlib.pyplot as plt

from .pairs import PairsFile, shift_line_median
from .plotting import parse_region, get_epigenetics, ContactMap, apply_matplotlib_style, format_ticks
from .neighbor import normalize_contact_map_neighbor
from . import gtf as gtf_loader
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def open_text(path: str):
    if path == '-':
        return sys.stdin
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def build_length_position_histogram(pairs_fp,
                                    chrom: str,
                                    start: int,
                                    end: int,
                                    length_min: int,
                                    length_max: int,
                                    bin_size: int,
                                    pair_type_filter: Optional[str] = 'UU') -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Stream a pairs file and accumulate a 2D histogram of counts over (x position bin, alignment length).

    For each read end, mark the full aligned span between its 5' and 3' mapped positions
    (pos5x..pos3x). Each covered x bin within the span increments the count at the
    corresponding alignment length row (|pos5 - pos3|).

    Returns (heatmap, x_bin_edges, y_bin_edges)
    """
    region_len = end - start
    num_x_bins = int(np.ceil(region_len / bin_size))
    # y bins are per bp from length_min..length_max inclusive
    y_values = np.arange(length_min, length_max + 1)
    heatmap = np.zeros((len(y_values), num_x_bins), dtype=np.int64)

    def _parse(value):
        return None if value == '.' else int(value)

    pf = PairsFile(pairs_fp)
    for line in pf:
        if pair_type_filter and line.get('pair_type', None) != pair_type_filter:
            continue

        # Process each end independently; use full alignment span [min(pos5,pos3), max(pos5,pos3)]
        for end_idx in (1, 2):
            chrom_key = f'chrom{end_idx}'
            pos5_key = f'pos5{end_idx}'
            pos3_key = f'pos3{end_idx}'

            if chrom_key not in line or line[chrom_key] != chrom:
                continue
            if pos5_key not in line or pos3_key not in line:
                continue

            p5 = _parse(line[pos5_key])
            p3 = _parse(line[pos3_key])
            if p5 is None or p3 is None:
                continue

            length_val = abs(p5 - p3)
            if length_val < length_min or length_val > length_max:
                continue

            span_start = max(min(p5, p3), start)
            span_end = min(max(p5, p3), end - 1)
            if span_end < span_start:
                continue

            s_bin = (span_start - start) // bin_size
            e_bin = (span_end - start) // bin_size
            s_bin = max(0, int(s_bin))
            e_bin = min(num_x_bins - 1, int(e_bin))
            if s_bin > e_bin:
                continue

            y_bin = length_val - length_min
            heatmap[y_bin, s_bin:e_bin + 1] += 1

    x_edges = np.arange(start, start + (num_x_bins + 1) * bin_size, bin_size)
    y_edges = np.arange(length_min, length_max + 2)  # +1 for inclusive edge
    return heatmap, x_edges, y_edges


def find_bigwigs(tracks_dir: Optional[str]) -> Dict[str, str]:
    if not tracks_dir:
        return {}
    paths = []
    paths += glob.glob(os.path.join(tracks_dir, '*.bw'))
    paths += glob.glob(os.path.join(tracks_dir, '*.bigWig'))
    tracks = {}
    for p in sorted(paths):
        name = os.path.basename(p)
        tracks[name] = p
    return tracks


@click.command()
@click.argument('region')
@click.option('--pairs', multiple=True, required=True, help='Input pairs file(s) (.pairs or .pairs.gz). Provide multiple --pairs to merge replicates.')
@click.option('--output', required=True, help='Output figure path (.svg recommended; .png/.pdf also supported).')
@click.option('--length-min', default=20, show_default=True, help='Minimum alignment length (bp), inclusive.')
@click.option('--length-max', default=80, show_default=True, help='Maximum alignment length (bp), inclusive.')
@click.option('--bin-size', default=1, show_default=True, help='Heatmap X-axis bin size (bp). Use 1 for single-bp resolution.')
@click.option('--pair-type', default='UU', show_default=True, help='Filter to this pair_type; empty to disable.')
@click.option('--tracks-dir', default=None, help='Directory containing bigWig tracks to plot as annotations.')
@click.option('--npz', default=None, help='Region npz file with balanced contact map for density (preferred).')
@click.option('--condition-name', default=None, help='Condition label for coloring (anatelo, midG1, prometa).')
@click.option('--title', default=None, help='Optional figure title.')
@click.option('--gtf', default=None, help='Optional GTF/GFF file to render a gene track.')
@click.option('--genes', default=None, help='Optional comma-separated gene names to highlight (others dim).')
@click.option('--log-counts/--no-log-counts', default=False, show_default=True, help='Use log1p(counts) for heatmap color scale.')
@click.option('--colorbar/--no-colorbar', default=False, show_default=True, help='Show a compact colorbar that does not shift axes.')
def main(region, pairs, output, length_min, length_max, bin_size, pair_type, tracks_dir, npz, condition_name, title, gtf, genes, log_counts, colorbar):
    """
    Generate a length (y: bp) vs genomic position (x) heatmap for a REGION using a pairs file.

    The y-axis spans [length-min, length-max]. For each end in each provided pairs file, we mark the full
    aligned span between its 5' and 3' mapped positions, and accumulate counts across all files.

    The figure also includes:
    - Contact density (marginal) for REGION if a cooler path is provided
    - All bigWig tracks found in --tracks-dir for REGION
    """
    apply_matplotlib_style()

    chrom, start, end = parse_region(region)

    # Build heatmap across one or more pairs files
    total_heatmap = None
    x_edges = None
    y_edges = None
    for pairs_path in pairs:
        with open_text(pairs_path) as fp:
            hm, xe, ye = build_length_position_histogram(
                pairs_fp=fp,
                chrom=chrom,
                start=start,
                end=end,
                length_min=length_min,
                length_max=length_max,
                bin_size=bin_size,
                pair_type_filter=pair_type if pair_type else None,
            )
            if total_heatmap is None:
                total_heatmap = hm
                x_edges = xe
                y_edges = ye
            else:
                total_heatmap += hm
    heatmap = total_heatmap if total_heatmap is not None else np.zeros((length_max - length_min + 1, int(np.ceil((end-start)/bin_size))), dtype=np.int64)

    # Prepare tracks and (optional) contact density
    tracks = find_bigwigs(tracks_dir)

    contact_x = None
    contact_density = None
    contact_color_map = {
        'anatelo': '#1f77b4',
        'midG1': '#ff7f0e',
        'prometa': '#2ca02c',
    }
    color = contact_color_map.get(condition_name, 'black') if condition_name else 'black'

    if npz:
        # Load balanced contact map from npz and apply neighbor normalization for consistency
        cm = ContactMap.from_npz(npz)
        balanced = normalize_contact_map_neighbor(cm.contact_map.copy(), max_prob=10.0, neighbor_prob=1.0)
        cm = ContactMap(balanced, cm.chrom, cm.start, cm.end, cm.resolution)
        contact_x = cm.x()
        contact_density = cm.get_marginal()

    # Figure layout: contact density (if any), heatmap, tracks, optional gene track
    n_tracks = len(tracks)
    height_ratios = []
    axes_count = 0
    if contact_density is not None:
        height_ratios.append(0.8)
        axes_count += 1
    height_ratios.append(2.5)  # heatmap
    axes_count += 1
    height_ratios.extend([0.5] * n_tracks)
    if gtf:
        height_ratios.append(0.6)  # gene track
        axes_count += 1
    axes_count += n_tracks

    fig, axs = plt.subplots(axes_count, 1, figsize=(12, sum(height_ratios) + 0.5), sharex=True,
                            gridspec_kw={'height_ratios': height_ratios})
    ax_idx = 0

    # Contact density
    if contact_density is not None:
        ax = axs[ax_idx]
        ax.plot(contact_x, contact_density, c=color, lw=1.5)
        ax.set_ylabel('Contact\ndensity', rotation=0, ha='right', va='center')
        format_ticks(ax, y=False)
        ax_idx += 1

    # Heatmap panel
    ax_hm = axs[ax_idx]
    extent = (x_edges[0], x_edges[-1], y_edges[0], y_edges[-1])
    display = np.log1p(heatmap) if log_counts else heatmap
    im = ax_hm.imshow(display, aspect='auto', origin='lower', extent=extent, cmap='coolwarm', interpolation='nearest')
    if colorbar:
        cax = inset_axes(ax_hm, width='40%', height='5%', loc='upper right', borderpad=0.8)
        cb = fig.colorbar(im, cax=cax, orientation='horizontal')
        cb.set_label('log1p(counts)' if log_counts else 'counts')
    ax_hm.set_ylabel('Alignment length (bp)')
    ax_hm.set_xlim(start, end)
    format_ticks(ax_hm, y=True)
    ax_idx += 1

    # BigWig tracks (binned to 200bp inside get_epigenetics to match heatmap default)
    for name, path in tracks.items():
        ax = axs[ax_idx]
        x_vals, vals = get_epigenetics(path, chrom, start, end, smoothing=0, bin=(bin_size != 1))
        ax.fill_between(x_vals, np.zeros_like(vals), vals, color='green')
        ax.set_ylabel(name, rotation=0, ha='right', va='center')
        ax.get_yticklabels()[:1]
        format_ticks(ax, y=False)
        ax_idx += 1

    # Optional gene track
    if gtf:
        ax = axs[ax_idx]
        gene_df = gtf_loader.dataframe(gtf)
        # Filter to region
        mask = (gene_df['seqname'] == chrom) & (gene_df['start'] <= end) & (gene_df['end'] >= start)
        gene_df = gene_df.loc[mask]
        # Prefer 'gene' features; fall back to 'transcript' if needed
        if 'feature' in gene_df.columns:
            pref = gene_df.loc[gene_df['feature'] == 'gene']
            if len(pref):
                gene_df = pref
            else:
                gene_df = gene_df.loc[gene_df['feature'] == 'transcript']
        # Parse highlight list
        highlight = set()
        if genes:
            highlight = set([g.strip() for g in genes.split(',') if g.strip()])
        # Draw genes
        for _, row in gene_df.iterrows():
            g_start = max(int(row['start']), start)
            g_end = min(int(row['end']), end)
            if g_start >= g_end:
                continue
            name = row.get('gene_name') or row.get('transcript_id') or row.get('gene_id') or ''
            is_highlight = (name in highlight) if name else False
            ax.hlines(0.5, g_start, g_end, colors='black' if is_highlight else 'gray', lw=2 if is_highlight else 1)
            if name:
                ax.text((g_start + g_end) / 2, 0.6, name, ha='center', va='bottom', fontsize=7,
                        color='black' if is_highlight else 'gray')
        ax.set_ylim(0, 1)
        ax.set_ylabel('Genes', rotation=0, ha='right', va='center')
        format_ticks(ax, y=False)

    axs[-1].set_xlim(start, end)
    format_ticks(axs[-1], y=False)

    if title:
        fig.suptitle(title, y=0.995)

    fig.savefig(output)
    plt.close(fig)


if __name__ == '__main__':
    main()


