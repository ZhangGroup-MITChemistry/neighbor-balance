import numpy as np
import pandas as pd
import cooler
import gzip
import matplotlib.pyplot as plt

try:
    from cooler.reduce import merge_breakpoints
except ImportError:
    from cooler._reduce import merge_breakpoints


class PairsFile:
    """Read and write pairs files."""
    def __init__(self, fp):
        self.fp = fp
        self.header = self.read_header()
        self.columns = self.get_columns()

    def read_header(self):
        header = []
        for line in self.fp:
            header.append(line)
            if line.startswith('#columns:'):
                return header
            if not line.startswith('#'):
                raise ValueError('No columns annotation found.')
        raise ValueError('No data found.')

    def get_columns(self):
        assert self.header[-1].startswith('#columns:')
        return self.header[-1].strip().split()[1:]

    def __iter__(self):
        for line in self.fp:
            yield self.process_line(line)

    def process_line(self, line):
        line = line.split()
        return {k: v for k, v in zip(self.columns, line)}

    def format_line(self, entries):
        return '\t'.join(str(entries[k]) for k in self.columns)

    @classmethod
    def as_pandas(cls, path, **kwargs):
        if path.endswith('.gz'):
            o = gzip.open
        else:
            o = open
        with o(path, 'rt') as f:
            columns = cls(f).columns
        return pd.read_csv(path, comment='#', sep='\t', header=None, names=columns, **kwargs)


def shift_line(line, protected_over_2):
    if line['pos1'] != '.':
        line['pos1'] = int(line['pos1'])
        if line['strand1'] == '+':
            line['pos1'] += protected_over_2
        else:
            line['pos1'] -= protected_over_2
    if line['pos2'] != '.':
        line['pos2'] = int(line['pos2'])
        if line['strand2'] == '+':
            line['pos2'] += protected_over_2
        else:
            line['pos2'] -= protected_over_2
    return line


def line_is_valid(line, chrom_sizes, *, regions=None, direction_to_keep=None, min_distance=0):
    chrom1 = line['chrom1']
    chrom2 = line['chrom2']
    pos1 = int(line['pos1'])
    pos2 = int(line['pos2'])
    strand1 = line['strand1']
    strand2 = line['strand2']
    pair_type = line['pair_type']

    # Only retain uniquely mapping reads.
    if pair_type != 'UU':
        return False

    # Determine the direction of the read and filter if specified.
    if strand1 == '+' and strand2 == '-':
        direction = 'inward'
    elif strand1 == '-' and strand2 == '+':
        direction = 'outward'
    elif strand1 == '+' and strand2 == '+':
        direction = 'tandemexit'
    elif strand1 == '-' and strand2 == '-':
        direction = 'tandementry'
    else:
        raise ValueError(f'Invalid strands: {strand1}, {strand2}')

    if direction_to_keep is not None and direction != direction_to_keep:
        return False

    # Filter reads that are too close or outside the chromosome.
    if pos1 <= 0 or pos2 <= 0 or pos1 > chrom_sizes[chrom1] or pos2 > chrom_sizes[chrom2]:
        return False
    if pos2 - pos1 < min_distance:
        return False

    # Filter reads that are outside the specified regions.
    if regions is not None:
        if chrom1 not in regions or chrom2 not in regions:
            return False
        if not regions[chrom1].overlaps(pos1) or not regions[chrom2].overlaps(pos2):
            return False
    return True


def combine_pixel_iterator(full, inward, bins, k):
    indexes = [c.open("r")["indexes/bin1_offset"] for c in [full, inward]]
    bin1_partition, cum_nrecords = merge_breakpoints(indexes, 1_000_000)
    starts = [0, 0]
    for bin1_id in bin1_partition[1:]:
        stops = [index[bin1_id] for index in indexes]

        pixels = full.pixels()[starts[0]:stops[0]].set_index(['bin1_id', 'bin2_id'])
        inward_pixels = inward.pixels()[starts[1]:stops[1]].set_index(['bin1_id', 'bin2_id'])
        pixels = pixels.join(inward_pixels, how='left', rsuffix='_inward')
        pixels['count_inward'] = pixels['count_inward'].fillna(0)

        chrom1 = bins.loc[pixels.index.get_level_values('bin1_id'), 'chrom'].to_numpy()
        chrom2 = bins.loc[pixels.index.get_level_values('bin2_id'), 'chrom'].to_numpy()
        start1 = bins.loc[pixels.index.get_level_values('bin1_id'), 'start'].to_numpy()
        start2 = bins.loc[pixels.index.get_level_values('bin2_id'), 'start'].to_numpy()

        diag_mask = (chrom1 == chrom2) & (np.abs(start2 - start1) <= k * full.binsize)

        pixels.loc[diag_mask, 'count'] -= pixels.loc[diag_mask, 'count_inward']
        pixels.loc[diag_mask, 'count'] = (4 * pixels.loc[diag_mask, 'count']) // 3
        pixels = pixels.reset_index().loc[:, ['bin1_id', 'bin2_id', 'count']]

        yield {k: v.values for k, v in pixels.items()}
        starts = stops


def remove_inward_reads_from_cooler(full_cool, inward_cool, output_cool, k):
    assert full_cool.endswith('.cool')
    assert inward_cool.endswith('.cool')
    assert output_cool.endswith('.cool')
    assert k >= 0

    full = cooler.Cooler(full_cool)
    bins = full.bins()[::]

    inward = cooler.Cooler(inward_cool)
    inward_bins = inward.bins()[::]

    assert bins['chrom'].cat.categories.equals(inward_bins['chrom'].cat.categories)
    assert np.all(bins['chrom'] == inward_bins['chrom'])
    assert np.all(bins['start'] == inward_bins['start'])
    assert np.all(bins['end'] == inward_bins['end'])

    assert full.binsize is not None
    assert full.binsize == inward.binsize

    cooler.create_cooler(output_cool, bins, combine_pixel_iterator(full, inward, bins, k))


####################################################################################################

def sliding_mean(data_array, window):
    new_list = []
    for i in range(data_array.shape[0]):
        if i <= window:
            indices = range(0, 2*i+1)
        else:
            indices = range(max(i - window, 0),
                            min(i + window + 1, len(data_array)))
        avg = 0
        for j in indices:
            avg += data_array[j]
        avg /= float(len(indices))
        new_list.append(avg)
    return np.array(new_list)


def count(data, high):
    data = data[(0 <= data) & (data < high)].astype(int)
    y = np.zeros(high)
    np.add.at(y, data, 1)
    return y


def ligation_distance_plot(counts, fname, colors, window=5):
    f, axs = plt.subplots(1, 2, figsize=(15, 5))
    for direction, _counts in counts.items():
        for ax in axs:
            ax.scatter(range(len(_counts)), _counts, s=5, alpha=0.5, c=colors[direction])
            ax.plot(range(len(_counts)), sliding_mean(_counts, window), c=colors[direction], label=direction)

    axs[1].set_yscale('log')
    axs[1].legend()
    axs[0].set_ylabel('read count')
    axs[1].set_ylabel('log read count')
    axs[0].set_xlabel('genomic separation (bp)')
    axs[1].set_xlabel('genomic separation (bp)')
    plt.savefig(fname)
    plt.close()


def windowed_distance_plot(df, window=200, high=10, pos1='pos1', pos2='pos2', ax=None, label=None):
    counts = np.zeros(high)
    for direction in df.direction.unique():
        _df = df.loc[df.direction == direction]
        dist = (_df[pos2]//window - _df[pos1]//window).to_numpy()
        counts += count(dist, high)
    counts /= counts[1]
    if ax is None:
        f, ax = plt.subplots()

    ax.scatter(range(high), counts, s=5, alpha=0.5)
    ax.plot(range(high), counts, label=label)
    ax.legend()
    ax.set_ylabel('relative read frequency')
    ax.set_xlabel(f'genomic separation (x {window} bp)')
    ax.set_ylim(0, 1)


def all_pairs_plots(pairs_fname, protected_over_2, nrows):
    colors = {'inward': 'cornflowerblue',
              'outward': 'mediumvioletred',
              'tandem_entry': 'darkorange',
              'tandem_exit': 'gold'}

    high = 1_000
    counts = {k: np.zeros(high) for k in ['inward', 'outward', 'tandem_entry', 'tandem_exit']}
    counts_shifted = {k: np.zeros(high) for k in ['inward', 'outward', 'tandem_entry', 'tandem_exit']}
    for df in PairsFile.as_pandas(pairs_fname, nrows=nrows, chunksize=1_000_000):
        df = df.loc[df.pair_type == 'UU']
        intra_chrom = df.chrom1 == df.chrom2
        print(f'frac intrachromosome  {sum(intra_chrom) / len(intra_chrom)}')
        df = df.loc[intra_chrom]
        near = (df.pos2 - df.pos1) < 1e5
        print(f'frac near  {sum(near) / len(near)}')

        df['direction'] = ''
        df.loc[(df.strand1 == '+') & (df.strand2 == '-'), 'direction'] = 'inward'
        df.loc[(df.strand1 == '-') & (df.strand2 == '+'), 'direction'] = 'outward'
        df.loc[(df.strand1 == '+') & (df.strand2 == '+'), 'direction'] = 'tandem_entry'
        df.loc[(df.strand1 == '-') & (df.strand2 == '-'), 'direction'] = 'tandem_exit'

        df['pos1_shifted'] = df.pos1
        df['pos2_shifted'] = df.pos2
        df.loc[df.direction == 'inward', 'pos1_shifted'] += protected_over_2
        df.loc[df.direction == 'inward', 'pos2_shifted'] -= protected_over_2
        df.loc[df.direction == 'outward', 'pos1_shifted'] -= protected_over_2
        df.loc[df.direction == 'outward', 'pos2_shifted'] += protected_over_2

        for direction in counts:
            _df = df.loc[df.direction == direction]
            dist = (_df.pos2 - _df.pos1).to_numpy()
            dist = dist[dist > 0]
            counts[direction] += count(dist, high)

        for direction in counts:
            _df = df.loc[df.direction == direction]
            dist = (_df.pos2_shifted - _df.pos1_shifted).to_numpy()
            dist = dist[dist > 0]
            counts_shifted[direction] += count(dist, high)

        # thresh = 1000
        # close = df.loc[(df.pos2 - df.pos1) <= thresh]
        # close_counts = [np.mean(close.direction == direction) for direction in colors]
        # far = df.loc[(df.pos2 - df.pos1) > thresh]
        # far_counts = [np.mean(far.direction == direction) for direction in colors]

    # f, ax = plt.subplots(1, 2, figsize=(10, 5))
    # ax[0].pie(close_counts, labels=colors.keys(), colors=colors.values(), normalize=True, autopct='%1.1f%%')
    # ax[0].set_title('genomic separation <= 1kb')
    # ax[1].pie(far_counts, labels=colors.keys(), colors=colors.values(), normalize=True, autopct='%1.1f%%')
    # ax[1].set_title('genomic separation > 1kb')
    # plt.savefig(f'{pairs_fname}_pie.png')
    # plt.close()

    ligation_distance_plot(counts, f'{pairs_fname}_raw.png', colors)
    ligation_distance_plot(counts_shifted, f'{pairs_fname}_shifted.png', colors)

    # f, ax = plt.subplots()
    # windowed_distance_plot(df, ax=ax, label='raw')
    # windowed_distance_plot(df, ax=ax, label='shifted', pos1='pos1_shifted', pos2='pos2_shifted')
    # plt.savefig(f'{pairs_fname}_windowed.png')
    # plt.close()
