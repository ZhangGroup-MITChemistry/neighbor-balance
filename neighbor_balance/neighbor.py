import numpy as np
from scipy.ndimage import gaussian_filter1d
import cooler
import logging


def get_neighbor_factors(diag, eps=0.0, average='arithmetic', main_diag=None):
    """
    Compute the neighbor factors for each window in a contact map.

    Parameters
    ----------
    diag: np.ndarray
        The first diagonal of the contact map: np.diagonal(contact_map, 1).
    eps: float
        The minimum value for the neighbor factors as a fraction of the mean neighbor factor.
        This is used to prevent pathological behavior when the neighbor factors are very small.
    average: str
        The method to average the neighbor factors. Options are 'harmonic', 'geometric', 'arithmetic', and 'mse'.
        The 'arithmetic' average is recommended.

    Returns
    -------
    np.ndarray
        The neighbor factors.
    """
    neighbor_factors = np.zeros(len(diag) + 1)

    if eps > 0:
        scaled_eps = np.nanmean(diag) * eps
        scaled_eps_inverse = np.nanmean(diag) / eps
        diag[diag < scaled_eps] = scaled_eps
        diag[diag > scaled_eps_inverse] = scaled_eps_inverse

    for i in range(len(diag) + 1):
        if i == 0:
            idx = [0]
        elif i == len(diag):
            idx = [-1]
        else:
            idx = [i-1, i]

        neighbors = diag[idx]
        if main_diag is not None:
            assert average == 'arithmetic', 'main_diag can only be used with arithmetic average'
            neighbors = np.append(neighbors, main_diag[i])

        if np.all(np.isnan(neighbors)):
            neighbor_factors[i] = np.nan
        else:
            if average == 'harmonic':
                neighbor_factors[i] = 1/np.nanmean(1/neighbors)
            elif average == 'geometric':
                prod = np.nanprod(neighbors)
                if not np.any(np.isnan(neighbors)):
                    prod = prod ** (1/len(neighbors))
                neighbor_factors[i] = prod
            elif average == 'arithmetic':
                neighbor_factors[i] = np.nanmean(neighbors)
            elif average == 'mse':
                neighbor_factors[i] = np.nansum(neighbors**2) / np.nansum(neighbors)
            else:
                raise ValueError(f'Unknown average method: {average}')
    return neighbor_factors


def normalize_contact_map_average(contact_map, max_prob=10.0, neighbor_prob=1.0):
    """
    Normalize the contact map by dividing by the average of the diagonal (a single scaler).

    This has no effect on the overall structure of the contacts, but just scales them.

    Parameters
    ----------
    contact_map: np.ndarray
        The contact map to be normalized.
    max_prob: float
        The maximum probability for a contact.
    neighbor_prob: float
        The probability of a contact between neighboring loci.

    Returns
    -------
    np.ndarray
        The normalized contact map.
    """
    average_neighbor = np.nanmean(np.diagonal(contact_map, 1))
    correction = neighbor_prob / average_neighbor
    contact_map = contact_map * correction
    mask = contact_map > max_prob
    count = np.sum(np.triu(mask, k=2))
    logging.info(f'Capping {count} non-diagonal contact probabilities to {max_prob}.')
    contact_map[mask] = max_prob
    return contact_map


def normalize_contact_map_neighbor(contact_map, bw=1, max_prob=10.0, neighbor_prob=1.0, max_iter=1, tol=1e-6, eps=1e-1,
                                   average='arithmetic', use_main_diag=False):
    """
    Balance the contact map such that the first off-diagonal elements are close to 1.

    For each locus, we compute a "neighbor factor" which is the average of the reciprocal of the contact probabilities
    with the neighboring loci. Each contact probability is then divided by the geometric mean of the neighbor factors
    for the two loci involved.

    Parameters
    ----------
    contact_map: np.ndarray
        The contact map to be normalized.
    bw: int
        The bandwidth of the Gaussian filter applied to the neighbor factors. Smoothing can help lower the variance if
        the first off-diagonal elements are noisy.
    max_prob: float
        The maximum probability for a contact.
    neighbor_prob: float
        The probability of a contact between neighboring loci.
    max_iter: int
        If greater than 1, iteratively perform neighbor balancing. Empirically, this seems to make the contact map
        very high variance and is not recommended.
    tol: float
        The tolerance for convergence if max_iter > 1.
    eps: float
        The minimum value for the neighbor factors as a fraction of the mean neighbor factor
         This is used to prevent pathological behavior when the neighbor factors are very small.
    average: str
        The method to average the neighbor factors. Options are 'harmonic', 'geometric', 'arithmetic', and 'mse'.
        The 'arithmetic' average is recommended.

    Returns
    -------
    np.ndarray
        The normalized contact map.
    """
    for i in range(max_iter):
        neighbor_factors = get_neighbor_factors(np.diagonal(contact_map, 1).copy(), eps=eps, average=average,
                                                main_diag=np.diagonal(contact_map, 0) if use_main_diag else None)
        if bw > 0:
            neighbor_factors = gaussian_filter1d(neighbor_factors, bw, mode='reflect')
        norm = np.sqrt(neighbor_factors.reshape(1, -1) * neighbor_factors.reshape(-1, 1))
        contact_map = contact_map / norm

        logging.info(f'Iteration: {i}, Mean: {np.mean(neighbor_factors):.2e}, Stddev: {np.std(neighbor_factors):.2e}')
        if np.std(neighbor_factors) < tol:
            break
    else:
        logging.warning('Neighbor balancing did not converge. This is expected.')
    return normalize_contact_map_average(contact_map, max_prob=max_prob, neighbor_prob=neighbor_prob)


def get_diagonal_for_chrom(clr, chrom, batch_size=1_000_000, balance='weight'):
    """
    Get the first diagonal of the contact map for the given chromosome from a cooler file.

    Parameters
    ----------
    clr : cooler.Cooler
        Cooler object.
    chrom : str
        Chromosome name.
    batch_size : int
        Batch size to calculate the neighbors.

    Returns
    -------
    np.ndarray
        Array with the average of the neighbors for each bin in the chromosome.
    """
    # I'm sure there is a better way to do this...
    assert batch_size <= 1_000_000, 'batch size must be less than 1 Mb.'
    batch_size = (batch_size // clr.binsize) * clr.binsize  # Make sure batch size is a multiple of the binsize.

    start, end = 0, clr.chromsizes[chrom]  # Half open, 0-based interval.

    # When a region is fetched, if any part of the region overlaps a bin, the whole bin is returned.
    # Therefore, for the whole chromosome, we will get ceil(end / clr.binsize) bins. The first diagonal
    # will then have one fewer elements.
    diagonal = np.zeros(np.ceil(end / clr.binsize).astype(int) - 1)

    for i, s in enumerate(range(start, end, batch_size)):
        # s and e are 0-based, half-open coordinates and should always be multiples of the binsize
        e = min(end, s + batch_size + clr.binsize)  # Add binsize to get overlap between batches.

        # These are the indices in the diagonal that will be written in this batch.
        s_i = s // clr.binsize
        e_i = np.ceil(e / clr.binsize).astype(int) - 1

        if not i % 10:
            logging.info(f'Getting neighbors for {chrom}:{s}')

        # Region coordinates are inclusive, 1-based, so we should add one to the start.
        cmap = clr.matrix(balance=balance).fetch(f'{chrom}:{s+1}-{e}')  
        diagonal[s_i:e_i] = np.diagonal(cmap, 1)
    return diagonal


def add_neighbor_factors_to_cooler(cool_fname, neighbor_res=200, batch_size=1_000_000, eps=10, overwrite=False):
    """
    Add neighbor balanced weights to the cooler file for each resolution.

    The new weights are stored in the 'weight_neighbor' field of the cooler file.
    """
    WEIGHT_NEIGHBOR = 'weight_neighbor'

    # Check if the cooler file already has a neighbor weight. if `overwrite`, remove it, otherwise raise an error.
    for res_path in cooler.fileops.list_coolers(cool_fname):
        clr = cooler.Cooler(f'{cool_fname}::{res_path}')
        if WEIGHT_NEIGHBOR in clr.bins().keys():
            if overwrite:
                logging.info(f'Overwriting {WEIGHT_NEIGHBOR} in {res_path}.')
                with clr.open("r+") as grp:
                    del grp["bins"][WEIGHT_NEIGHBOR]
            else:
                raise ValueError(f'Cooler {res_path} already has a neighbor weight. Use --overwrite to overwrite.')

    # Compute the neighbor factors for the neighbor resolution.
    clr = cooler.Cooler(f'{cool_fname}::/resolutions/{neighbor_res}')
    all_neighbors = {}
    for chrom in clr.chromsizes.index:
        diagonal = get_diagonal_for_chrom(clr, chrom, batch_size=batch_size)
        neighbors = get_neighbor_factors(diagonal)
        assert len(neighbors) == np.ceil(clr.chromsizes[chrom] / neighbor_res)
        all_neighbors[chrom] = neighbors

    # Compute the neighbor factors for each resolution, as the average over the bins, and add them to the cooler file.
    for res_path in cooler.fileops.list_coolers(cool_fname):
        clr = cooler.Cooler(f'{cool_fname}::{res_path}')
        assert clr.binsize % neighbor_res == 0, 'Resolution must be a multiple of the neighbor resolution.'
        cg_neighbors = []
        for chrom in clr.chromsizes.index:
            bin_count = clr.binsize // neighbor_res
            for i in range(0, len(all_neighbors[chrom]), bin_count):
                vals = all_neighbors[chrom][i:i + bin_count]
                if np.all(np.isnan(vals)) or len(vals) == 0:
                    # Use the whole chromosome mean if none of the neighbors in the window are valid.
                    # It should be rare that the row/column is unmasked and the neighbors are all nan,
                    # but it can happen. Using the mean of the whole chromosome is a reasonable default.
                    cg_neighbors += [np.nanmean(all_neighbors[chrom])]
                else:
                    cg_neighbors += [np.nanmean(vals)]
        cg_neighbors = np.array(cg_neighbors)

        # Apply eps to the cg_neighbors to prevent pathological behavior.
        median = np.nanmedian(cg_neighbors)
        cg_neighbors = np.clip(cg_neighbors, median / eps, median * eps)

        weights = clr.bins()['weight'][:] / np.sqrt(cg_neighbors)
        with clr.open("r+") as grp:
            grp["bins"].create_dataset(WEIGHT_NEIGHBOR, data=weights)
