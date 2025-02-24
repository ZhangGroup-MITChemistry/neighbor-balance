import pytest
from neighbor_balance.neighbor import *
import numpy as np
import pandas as pd


def make_synthetic_cooler(contact_map, name, noise=False):
    """
    Create a cooler file with the given contact map.
    """
    if noise:
        a = np.random.lognormal(0, 1, contact_map.shape[0])
        contact_map = contact_map * a.reshape(-1, 1) * a.reshape(1, -1)

    bins = pd.DataFrame([('chr1', i, i + 200) for i in range(0, 200 * contact_map.shape[0] + 1, 200)],
                        columns=['chrom', 'start', 'end'])

    n = 10_000_000
    pixels = []
    for i in range(contact_map.shape[0]):
        for j in range(i, contact_map.shape[0]):
            pixels.append((i, j, (n * contact_map[i, j]).astype(int)))
    pixels = pd.DataFrame(pixels, columns=['bin1_id', 'bin2_id', 'count'])

    cooler.create_cooler(f'{name}.cool', bins, pixels)
    cool = cooler.Cooler(f'{name}.cool')
    weights, info = cooler.balance_cooler(cool, ignore_diags=1, mad_max=0, min_nnz=1)
    bins['weight'] = weights
    cooler.create_cooler(f'{name}.mcool::/resolutions/200', bins, pixels)


def test_get_neighbor_factors():
    diag = np.array([1, 2, 3, 4])
    eps = 0
    average = 'harmonic'
    assert np.allclose(get_neighbor_factors(diag, eps, average), [1.0, 2/(1/1+1/2), 2/(1/2+1/3), 2/(1/3+1/4), 4.0])


def test_add_to_cooler():
    contact_map = np.array([[1, 2],
                            [2, 1]])
    make_synthetic_cooler(contact_map, 'test')

    diagonal = get_diagonal_for_chrom(cooler.Cooler('test.mcool::/resolutions/200'), 'chr1')
    assert np.allclose(diagonal, [1, np.nan], equal_nan=True)

    neighbors = get_neighbor_factors(diagonal)
    assert np.allclose(neighbors, [1, 1, np.nan], equal_nan=True)

    add_neighbor_factors_to_cooler('test.mcool')
    cool = cooler.Cooler('test.mcool::/resolutions/200')
    df = cool.bins()[:]
    assert np.allclose(df['weight_neighbors'], neighbors, equal_nan=True)
    assert np.allclose(df['weight_neighbors_times_ice'], df['weight'] * neighbors, equal_nan=True)


def test_neighbor_factor_zero():
    # Contact map that gives a neighbor factor of zero
    contact_map = [[0.016,  0.069,  0.0,    np.nan, 0.407],
                   [0.069,  0.0,    0.0,    np.nan, 0.405],
                   [0.0,    0.0,    0.0,    np.nan, 0.534],
                   [np.nan, np.nan, np.nan, np.nan, np.nan],
                   [0.407, 0.405,   0.534,  np.nan, 0.032]]
    contact_map = np.array(contact_map)

    nf = get_neighbor_factors(contact_map, eps=0)
    assert nf[2] == 0

    eps = 1e-1
    nf = get_neighbor_factors(contact_map, eps=eps)
    assert nf[2] > 0
