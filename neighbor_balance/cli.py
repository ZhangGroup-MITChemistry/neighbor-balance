import click
import logging
import sys
import pandas as pd
import numpy as np
from intervaltree import IntervalTree, Interval
from .pairs import PairsFile, shift_line, line_is_valid, remove_inward_reads_from_cooler, all_pairs_plots, get_base_ps
from .neighbor import add_neighbor_factors_to_cooler, normalize_contact_map_neighbor
from .plotting import ContactMap, parse_region
from .ice import ice_balance_with_interpolation, get_capture_rates


@click.group()
def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")


@main.command()
@click.option('--protected-over-2', default=65, help='Number of bases to protect around the nucleosome center.')
def shift_pairs(protected_over_2):
    """
    Shift the positions by protected-over-2 in the 3' direction of the read.
    """
    pf = PairsFile(sys.stdin)
    sys.stdout.write(''.join(pf.header))
    for line in pf:
        line = shift_line(line, protected_over_2)
        sys.stdout.write(pf.format_line(line) + '\n')


@main.command()
@click.argument('chrom_sizes')
@click.option('--regions', default=None, help='Bed file with regions to keep.')
@click.option('--direction', default=None, help='Only keep reads in this direction.')
@click.option('--min-distance', default=0, help='Filter reads with nucleosomes closer than this distance.')
def filter_pairs(chrom_sizes, regions, direction, min_distance):
    """
    Filter pairs file.

    The following filters are applied:

    If regions is specified, only keep reads where both ends are in the specified regions.

    If direction is specified, only keep reads in the specified direction.
        Must be one of ['inward', 'outward', 'tandementry', 'tandemexit'].
    """
    chrom_sizes = pd.read_csv(chrom_sizes, sep='\t', header=None, names=['chrom', 'size'])
    chrom_sizes = {k: v for k, v in zip(chrom_sizes['chrom'], chrom_sizes['size'])}

    if regions is not None:
        regions_df = pd.read_csv(regions, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'])
        regions = {}
        for i, row in regions_df.iterrows():
            if row['chrom'] not in regions:
                regions[row['chrom']] = IntervalTree()
            regions[row['chrom']].add(
                Interval(row['start'] + 1, row['end'] + 1))  # bed is 0-based, half-open. pairs are 1-based.

    pf = PairsFile(sys.stdin)
    sys.stdout.write(''.join(pf.header))
    for line in pf:
        if line_is_valid(line, chrom_sizes, regions=regions, direction_to_keep=direction, min_distance=min_distance):
            sys.stdout.write(pf.format_line(line) + '\n')


@main.command()
@click.argument('pairs_fname')
@click.option('--protected-over-2', default=65, help='Number of bases to protect around the nucleosome center.')
@click.option('--nrows', default=1_000_000, help='Number of reads to analyze.')
def analyze_pairs(pairs_fname, protected_over_2, nrows):
    """
    Analyze pairs file.
    """
    all_pairs_plots(pairs_fname, protected_over_2, nrows)


@main.command()
@click.argument('pairs_fname')
@click.option('--nrows', default=None, type=int, help='Number of reads to analyze.')
@click.option('--max-distance', default=100_000, help='Maximum distance to consider.')
@click.option('--output', default=None, help='Output file name.')
def base_ps(pairs_fname, nrows, max_distance, output):
    """
    Save base resolution P(s) curve to a csv file.
    """
    counts = get_base_ps(pairs_fname, nrows=nrows, high=max_distance)
    counts = pd.DataFrame(counts)
    counts['distance'] = np.arange(len(counts))

    if output is None:
        output = pairs_fname + '.base_ps.csv'
    counts.to_csv(output, index=False)


@main.command()
@click.argument('full_cool')
@click.argument('inward_cool')
@click.argument('output_cool')
@click.option('--k', default=2, help='Number of diagonals to subtract the inward cool from the full cool.')
def remove_inward(full_cool, inward_cool, output_cool, k):
    """
    Remove reads counts in INWARD_COOL from the counts in FULL_COOL for the first K diagonals
    and scale the remaining counts by 4/3 and write to OUTPUT_COOL.
    """
    remove_inward_reads_from_cooler(full_cool, inward_cool, output_cool, k)


@main.command()
@click.argument('cool_fname')
@click.option('--neighbor-res', default=200, help='Resolution of the contact map used for computing neighbor weights.'
              ' This should be between 200 and 400 bps. All other resolutions in the cooler should be multiples of this value.')
@click.option('--batch-size', default=1_000_000, help='Batch size to calculate the neighbors.')
def neighbor_balance_cooler(cool_fname, neighbor_res, batch_size):
    """
    Add the inverse of the neighbors to the cooler.
    """
    add_neighbor_factors_to_cooler(cool_fname, neighbor_res=neighbor_res, batch_size=batch_size)


@main.command()
@click.argument('output_fname')
@click.argument('cool-fname')
@click.argument('region')
@click.option('--resolution', default=200, help='Resolution of the contact map.')
@click.option('--capture-probes-path', default=None, help='Path to the capture probes bed file.')
@click.option('--mad-max', default=2, help='Minimum coverage in terms of medium absolute deviations.')
@click.option('--min-nnz', default=100, help='Minimum number of non-zero values for the contact map.')
@click.option('--min-pair-capture-rate', default=0.2, help='Minimum pair capture rate for each pixel.')
@click.option('--max-iter', default=50, help='Maximum number of iterations for the balancing algorithm.')
@click.option('--tol', default=1e-5, help='Tolerance for the balancing algorithm.')
@click.option('--correct-for-flanks', default=True, help='Correct for flanking regions.')
@click.option('--sigma-scale', default=2.0, help='Scale for the smoothing function.')
@click.option('--sigma-exponent', default=0.2, help='Exponent for the smoothing function.')
@click.option('--sigma-plateau', default=2, help='Plateau for the smoothing function.')
def region_balance(output_fname, cool_fname, region, resolution, capture_probes_path, **normalization_params):
    """
    Balance individual regions of the contact map and save as a numpy file.

    This routine is intended to be used with region-capture data. It provides balancing with interpolation of missing
    values, accounts for flanking regions, and masks values not covered by a capture probe on either side. This
    functionality is not currently supported for genome-wide cooler files, but it could be added in the future.

    Note that this function does not perform neighbor balancing! Neighbor balancing on small regions is very fast and
    there is no need to do it in advance.

    To load saved contact map and subsequently neighbor balance it, use the following code:
    >>> loaded = np.load('output_fname.npz', allow_pickle=True)
    >>> contact_map = loaded['array']
    >>> metadata = loaded['metadata'].item()
    >>> contact_map = normalize_contact_map_neighbor(contact_map)
    """
    config = {
        'contact_map_path': cool_fname,
        'region': region,
        'resolution': resolution,
        'capture_probes_path': capture_probes_path,
        'normalization_params': normalization_params
    }

    # Load the contact map. Have contacts point to the underlying numpy array; we don't need the extra functionality.
    contacts = ContactMap.from_cooler(path=config['contact_map_path'], resolution=config['resolution'],
                                      region=config['region'], balance=False)
    contacts = contacts.contact_map

    # Load capture rates.
    if 'capture_probes_path' in config and config['capture_probes_path'] is not None:
        chrom, start, end = parse_region(config['region'])
        capture_rates = get_capture_rates(config['capture_probes_path'],
                                          chrom, start, end, bin_size=config['resolution'])
    else:
        capture_rates = None

    # Balance and save.
    contacts = ice_balance_with_interpolation(contacts, capture_rates=capture_rates, **config['normalization_params'])
    np.savez(output_fname, array=contacts, metadata=config)
