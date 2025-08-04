import pytest
from neighbor_balance.pairs import remove_inward_reads_from_cooler, line_is_valid, shift_line, PairsFile
import numpy as np
import cooler
import pandas as pd


def write_cool_file(contact_map, filename, start=0, chrom='chr1', switch=1e10, start2=0, chrom2='chr1'):
    pixels = []
    for i in range(contact_map.shape[0]):
        for j in range(i, contact_map.shape[1]):
            pixels += [(i, j, contact_map[i, j])]

    pixels = pd.DataFrame(pixels, columns=['bin1_id', 'bin2_id', 'count'])
    bins = []
    for i in range(contact_map.shape[0]):
        if i < switch:
            bins += [(chrom, start+i*1000, start+(i+1)*1000)]
        else:
            bins += [(chrom2, start2+(i-switch)*1000, start2+(i-switch+1)*1000)]

    bins = pd.DataFrame(bins, columns=['chrom', 'start', 'end'])

    cooler.create_cooler(filename, bins, pixels)


def test_combine_subtract():
    full = np.ones((5, 5))
    inward = np.ones((5, 5))
    write_cool_file(full, 'full.cool')
    write_cool_file(inward, 'inward.cool')
    remove_inward_reads_from_cooler('full.cool', 'inward.cool', 'combined.cool', 2)
    mat = cooler.Cooler('combined.cool').matrix(balance=False, sparse=False)[:]

    assert np.diagonal(mat, 0).sum() == 0
    assert np.diagonal(mat, 1).sum() == 0
    assert np.diagonal(mat, 2).sum() == 0
    assert np.diagonal(mat, 3).sum() == 2
    assert np.diagonal(mat, 4).sum() == 1


def test_combine_scale():
    full = 4*np.ones((5, 5))
    inward = np.ones((5, 5))
    inward[0, 1] = 0
    write_cool_file(full, 'full.cool')
    write_cool_file(inward, 'inward.cool')
    remove_inward_reads_from_cooler('full.cool', 'inward.cool', 'combined.cool', 2)
    mat = cooler.Cooler('combined.cool').matrix(balance=False, sparse=False)[:]

    assert mat[0, 0] == int(3 * 4 / 3)  # main diagonal, subtract 1
    assert mat[0, 1] == int(4*4/3)  # first off-diagonal, subtract 0
    assert mat[1, 2] == int(3*4/3)  # second off-diagonal, subtract 1
    assert mat[0, 2] == int(3*4/3)  # first off-diagonal, subtract 1

    for i in range(full.shape[0]):
        for j in range(i+2, full.shape[1]):
            assert mat[i, j] == full[i, j]


def test_combine_mismatch_starts():
    full = np.ones((5, 5))
    inward = np.ones((5, 5))
    write_cool_file(full, 'full.cool')
    write_cool_file(inward, 'inward.cool', start=1)

    with pytest.raises(AssertionError):
        remove_inward_reads_from_cooler('full.cool', 'inward.cool', 'combined.cool', 2)


def test_combine_mismatch_chrom():
    full = np.ones((5, 5))
    inward = np.ones((5, 5))
    write_cool_file(full, 'full.cool')
    write_cool_file(inward, 'inward.cool', chrom='chr2')

    with pytest.raises(AssertionError):
        remove_inward_reads_from_cooler('full.cool', 'inward.cool', 'combined.cool', 2)


def test_combine_mismatch_switch():
    full = np.ones((5, 5))
    inward = np.ones((5, 5))
    write_cool_file(full, 'full.cool', switch=2)
    write_cool_file(inward, 'inward.cool', switch=3)

    with pytest.raises(AssertionError):
        remove_inward_reads_from_cooler('full.cool', 'inward.cool', 'combined.cool', 2)


def test_combine_disjoint():
    full = np.ones((5, 5))
    inward = np.ones((5, 5))
    write_cool_file(full, 'full.cool', switch=2, start2=1_000_000)
    write_cool_file(inward, 'inward.cool', switch=2, start2=1_000_000)
    remove_inward_reads_from_cooler('full.cool', 'inward.cool', 'combined.cool', 2)
    mat = cooler.Cooler('combined.cool').matrix(balance=False, sparse=False)[:]

    assert np.diagonal(mat, 0).sum() == 0
    assert np.diagonal(mat, 1).sum() == 1  # One pixel in the disjoint region
    assert np.diagonal(mat, 2).sum() == 2  # Two pixels in the disjoint region
    assert np.diagonal(mat, 3).sum() == 2  # Same as before
    assert np.diagonal(mat, 4).sum() == 1


def test_read_pairs():
    example_pairs = 'example.pairs'
    with open(example_pairs, 'r') as f:
        pf = PairsFile(f)
        for line in pf:
            pass


def make_line(read_id='.', chrom1='chr1', pos1=1000, chrom2='chr1', pos2=2000,
              strand1='+', strand2='+', read_type='UU', mapq1=10, mapq2=10):
    return {
        'read_id': read_id,
        'chrom1': chrom1,
        'pos1': pos1,
        'chrom2': chrom2,
        'pos2': pos2,
        'strand1': strand1,
        'strand2': strand2,
        'pair_type': read_type,
        'mapq1': mapq1,
        'mapq2': mapq2,
    }


def chrom_sizes():
    return {'chr1': 1_000_000, 'chr2': 1_000_000}


def test_line_processor():
    line = make_line()
    assert line_is_valid(line, chrom_sizes())


def test_line_processor_shift():
    pos1 = 1000
    pos2 = 2000
    strand1 = '+'
    strand2 = '-'

    line = make_line(pos1=pos1, pos2=pos2, strand1=strand1, strand2=strand2)
    new_line = shift_line(line, protected_over_2=70)
    assert new_line == make_line(pos1=pos1+70, pos2=pos2-70, strand1=strand1, strand2=strand2)

    strand1 = '+'
    strand2 = '+'
    line = make_line(pos1=pos1, pos2=pos2, strand1=strand1, strand2=strand2)
    new_line = shift_line(line, protected_over_2=70)
    assert new_line == make_line(pos1=pos1+70, pos2=pos2+70, strand1=strand1, strand2=strand2)

    strand1 = '-'
    strand2 = '+'
    line = make_line(pos1=pos1, pos2=pos2, strand1=strand1, strand2=strand2)
    new_line = shift_line(line, protected_over_2=70)
    assert new_line == make_line(pos1=pos1-70, pos2=pos2+70, strand1=strand1, strand2=strand2)

    strand1 = '-'
    strand2 = '-'
    line = make_line(pos1=pos1, pos2=pos2, strand1=strand1, strand2=strand2)
    new_line = shift_line(line, protected_over_2=70)
    assert new_line == make_line(pos1=pos1-70, pos2=pos2-70, strand1=strand1, strand2=strand2)


def test_line_processor_shift_dots():
    pos1 = 1000
    pos2 = 2000
    strand1 = '+'
    strand2 = '-'

    line = make_line(pos1='.', pos2=pos2, strand1=strand1, strand2=strand2)
    new_line = shift_line(line, protected_over_2=70)
    assert new_line == make_line(pos1='.', pos2=pos2-70, strand1=strand1, strand2=strand2)

    line = make_line(pos1=pos1, pos2='.', strand1=strand1, strand2=strand2)
    new_line = shift_line(line, protected_over_2=70)
    assert new_line == make_line(pos1=pos1+70, pos2='.', strand1=strand1, strand2=strand2)


def test_line_processor_filter_direction():
    line = make_line(strand1='+', strand2='-')
    assert line_is_valid(line, chrom_sizes(), direction_to_keep='inward')

    line = make_line(strand1='+', strand2='+')
    assert not line_is_valid(line, chrom_sizes(), direction_to_keep='inward')


def test_line_processor_filter_min_dist():
    line = make_line(strand1='+', strand2='-', pos1=1000, pos2=999)
    assert not line_is_valid(line, chrom_sizes())

def test_line_is_valid():
    # Intra-chromosomal contact with pos2 < pos1
    line = make_line(chrom1='chr1', pos1=2000, chrom2='chr2', pos2=1000, strand1='+', strand2='+')
    assert line_is_valid(line, chrom_sizes())

    # Intra-chromosomal contact with pos1 < pos2
    line = make_line(chrom1='chr1', pos1=1000, chrom2='chr2', pos2=2000, strand1='+', strand2='+')
    assert line_is_valid(line, chrom_sizes())
