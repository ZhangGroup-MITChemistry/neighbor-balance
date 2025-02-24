import click
import sys
import pandas as pd
from intervaltree import IntervalTree, Interval
from .pairs import PairsFile, shift_line, line_is_valid, remove_inward_reads_from_cooler, all_pairs_plots
from .neighbor import add_neighbor_factors_to_cooler


@click.group()
def main():
    pass


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
            sys.stdout.write(pf.format_line(line))


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
@click.option('--batch-size', default=1_000_000, help='Batch size to calculate the neighbors.')
def neighbor_balance_cooler(cool_fname, batch_size):
    """
    Add the inverse of the neighbors to the cooler.
    """
    add_neighbor_factors_to_cooler(cool_fname, batch_size=batch_size)
