import numpy as np
import pandas as pd
from .smoothing import interpolate_nans, smoothing_sigma, interpolate_diagonals, get_distance_average


def ice_balance_with_interpolation(contact_map,
                                   mad_max=5, min_nnz=10, min_pair_capture_rate=0.2, capture_rates=None,
                                   max_iter=50, tol=1e-5, correct_for_flanks=False,
                                   sigma_scale=2.0, sigma_exponent=0.2, sigma_plateau=2):
    """
    Balance the contact map with rows/entries filtered by various criteria and interpolation of missing values.

    `capture_rates` can be obtained using:
        `chrom, start, end = parse_region(region)`, then
        `get_capture_rates(capture_probes_bed, chrom, start, end, bin_size=bin_size)`.
    """
    assert contact_map.dtype == np.int32, \
        f'Contact map must be raw counts and thereby type int, not {contact_map.dtype}.'

    if capture_rates is None:
        print('No capture rates provided. Assuming uniform capture rates.', flush=True)
        capture_rates = np.ones(contact_map.shape[0])

    contact_map = filter_bins(contact_map.astype(float), capture_rates,
                              mad_max=mad_max, min_nnz=min_nnz, min_pair_capture_rate=min_pair_capture_rate)

    def sigma(k):
        return smoothing_sigma(k, scale=sigma_scale, exponent=sigma_exponent, plateau=sigma_plateau)

    contact_map, _ = ice_balance(contact_map, interpolate=False,
                                 max_iter=max_iter, tol=tol, correct_for_flanks=correct_for_flanks)
    contact_map, _ = ice_balance(contact_map, interpolate=True, sigma=sigma,
                                 max_iter=max_iter, tol=tol, correct_for_flanks=correct_for_flanks)
    assert ~np.any(np.isnan(contact_map)), 'ICE balancing did not interpolate missing values.'
    contact_map, _ = ice_balance(contact_map, interpolate=False,
                                 max_iter=max_iter, tol=tol, correct_for_flanks=correct_for_flanks)
    return contact_map


def get_capture_rates(probe_bed_file, chrom, start, end, bin_size=200):
    """
    Read a bed file with capture probe positions and compute the capture rates for each locus.

    Parameters
    ----------
    probe_bed_file: str
        The path to the bed file with the capture probes.
    chrom: str
        The chromosome of the region.
    start: int
        The start of the region.
    end: int
        The end of the region.
    bin_size: int
        The size of the bins for the capture rates.

    Returns
    -------
    np.ndarray
        The capture rates for each locus.
    """
    capture_probes = pd.read_csv(probe_bed_file, sep='\t', header=None, names=['chrom', 'start', 'end'])

    start_in_region = (start <= capture_probes['start']) & (capture_probes['start'] < end)
    end_in_region = (start <= capture_probes['end']) & (capture_probes['end'] < end)
    capture_probes = capture_probes[(capture_probes['chrom'] == chrom) & (start_in_region | end_in_region)]

    capture_rates = np.zeros((end - start) // bin_size)
    for _, row in capture_probes.iterrows():
        for pos in range(row['start'], row['end']):
            if start <= pos < end:
                capture_rates[(pos - start) // bin_size] += 1/bin_size
    return capture_rates


def get_marginal(x, k=1, correct_for_flanks=False):
    """
    Compute the marginals of the contact map.

    Parameters
    ----------
    x: np.ndarray
        The contact map.
    k: int
        Number of diagonals to exclude from the marginal.
    correct_for_flanks: bool
        Whether to add the average of the flanks to the marginal.

    Returns
    -------
    np.ndarray
        The marginals.
    """
    if k > 0:
        x = np.triu(x, k=k) + np.triu(x, k=k).T

    marginal = np.nansum(x, axis=0)

    if correct_for_flanks:
        n = len(marginal)
        ps = get_distance_average(x)
        ps = interpolate_nans(ps)
        for i in range(len(marginal)):
            if marginal[i] != 0:
                marginal[i] += ps[n-i:].sum()  # Add the average of the flanks to the right
                marginal[i] += ps[i:].sum()  # and to the left
    return marginal


def filter_bins(contact_map, capture_rates, mad_max=5, min_nnz=10, min_pair_capture_rate=0.2):
    """
    Filters rows/columns and entries in the contact map based on various criteria.

    Entries failing the filters are set to NaN.

    The default filters match the ones used in cooler (as of 08/2024).

    However, this implementation differs in that the filters are applied iteratively until no updates are made. This is
    necessary because there are cases where the capture rate filter puts a row/column below the other thresholds,
    which results in problems during ICE balancing.

    Parameters
    ----------
    contact_map: np.ndarray
        The contact map to be filtered.
    capture_rates: np.ndarray
        The capture rates for each locus.
    mad_max: float
        The maximum median absolute deviation for the MAD-max filter.
    min_nnz: int
        The minimum number of non-zero entries for the min-nnz filter.
    min_pair_capture_rate: float
        The minimum capture rate for a pair of loci.

    Returns
    -------
    np.ndarray
        The filtered contact map.
    """
    contact_map = contact_map.copy()

    # MAD-max filter.
    if mad_max > 0:
        log_marginal = np.log(get_marginal(contact_map))
        mad = np.median(np.abs(log_marginal - np.median(log_marginal)))
        mask = log_marginal < np.median(log_marginal) - mad_max * np.median(mad)
        all_nan = np.all(np.isnan(contact_map), axis=1)
        n_changes = np.sum(mask & ~all_nan)
        print(f'Removing {n_changes} rows/columns with low coverage by mad-max.')
        contact_map[mask] = np.nan

    # Pair capture rate filter
    if min_pair_capture_rate > 0:
        pair_capture_rates = (capture_rates.reshape(-1, 1) + capture_rates.reshape(1, -1))
        mask = pair_capture_rates < min_pair_capture_rate
        n_changes = np.sum(mask & ~np.isnan(contact_map))
        print(f'Removing {n_changes} entries with low capture rate.')
        contact_map[mask] = np.nan

    # Number of non-zero pixel filter.
    if min_nnz > 0:
        mask = np.sum(contact_map > 0, axis=1) < min_nnz
        all_nan = np.all(np.isnan(contact_map), axis=1)
        n_changes = np.sum(mask & ~all_nan)
        print(f'Removing {n_changes} rows/columns with low coverage by min-nnz.')
        contact_map[mask] = np.nan
        contact_map[:, mask] = np.nan

    return contact_map


def filtered_map_is_valid(contact_map):
    """
    Checks if the contact map is valid to be balanced.

    Specifically, the contact map must not have any marginals that are zero. If this occurs, it will result in a
    division by zero during ICE balancing.

    However, there can be masked rows/columns in the contact map that are all set to NaN. While this does technically
    result in a division by zero, the values being divided by zero are already NaN, so it is not a problem.

    Additionally, the contact map must be symmetric.
    """
    # Check if marginals are non-zero or all NaN.
    marginal = get_marginal(contact_map)
    nans = np.all(np.isnan(contact_map), axis=1)
    if not np.all((marginal > 0) | nans):
        return False

    # Check if NaNs are symmetric.
    nans = np.isnan(contact_map)
    if not np.all(nans == nans.T):
        return False

    # Check if non-NaNs are symmetric.
    if not np.allclose(contact_map[~nans], contact_map.T[~nans], atol=1e-4):
        return False

    return True


def ice_balance(contact_map, sigma=2.0, max_iter=50, tol=1e-5, interpolate=False, correct_for_flanks=False):
    """
    Ice balancing algorithm with nan interpolation.

    Important: it is expected that the contact map has been filtered before balancing.

    Parameters
    ----------
    contact_map: np.ndarray
        The contact map to be balanced.
    sigma: float or callable
        The standard deviation of the Gaussian filter used for interpolation.
    max_iter: int
        The maximum number of iterations.
    tol: float
        The tolerance for convergence of the variance.
    interpolate: bool
        Whether to interpolate NaN values or treat them as zero.
    correct_for_flanks: bool
        Whether to add the average of the flanks when computing the marginals.

    Returns
    -------
    np.ndarray
        The balanced contact map.
    """
    if not filtered_map_is_valid(contact_map):
        raise ValueError('Contact map is not valid for ICE balancing.')

    mask = np.isnan(contact_map)
    W = contact_map.copy()
    if interpolate:
        W = interpolate_diagonals(W, sigma=sigma)
    else:
        W[mask] = 0.0
    B = np.ones(W.shape[0])
    for i in range(max_iter):
        S = get_marginal(W, correct_for_flanks=correct_for_flanks)
        dB = S / np.mean(S[S != 0])
        dB[S == 0] = 1
        W /= dB.reshape(-1, 1) * dB.reshape(1, -1)
        W[mask] = np.nan
        if interpolate:
            W = interpolate_diagonals(W, sigma=sigma)
        B *= dB
        print(f'Iteration {i}: mean={np.mean(S[S !=0 ])}, variance={np.var(S[S !=0 ])}')
        if np.var(S) < tol:
            break
    else:
        print('Warning: ICE balancing did not converge.')

    if not interpolate:
        W[mask] = np.nan
    return W, B
