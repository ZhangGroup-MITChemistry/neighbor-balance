import numpy as np
from scipy.ndimage import gaussian_filter


def get_distance_average(contact_map, show_self=False):
    """
    Compute the average of the contact frequencies for each distance.
    """
    n = contact_map.shape[0]
    if not show_self:
        n -= 1
    y = np.zeros(n)
    for i in range(n):
        offset = i if show_self else i + 1
        diag = contact_map.diagonal(offset=offset)
        if np.all(np.isnan(diag)):
            y[i] = np.nan
        else:
            y[i] = np.nanmean(diag)
    return y


def interpolate_nans(x, sigma=0):
    """
    Interpolate NaN values using linear interpolation.

    If `sigma` is greater than 0, the (non-NaN) values are first smoothed with a Gaussian filter before interpolation
    to reduce the effect of noise.

    Parameters
    ----------
    x: np.ndarray
        The array to be interpolated.
    sigma: float
        The standard deviation of the Gaussian filter.

    Returns
    -------
    np.ndarray
        The interpolated array.
    """
    x = x.copy()
    nans = np.isnan(x)
    if sigma > 0:
        x_smooth = nan_gaussian_filter(x, sigma=sigma, nan_policy='mask')
    else:
        x_smooth = x
    x[nans] = np.interp(np.arange(len(nans))[nans], np.arange(len(nans))[~nans], x_smooth[~nans])
    return x


def nan_gaussian_filter(x, sigma, nan_policy='mask'):
    """
    Apply a Gaussian filter to an array, ignoring NaN values.

    Parameters
    ----------
    x: np.ndarray
        The array to be smoothed.
    sigma: float
        The standard deviation of the Gaussian kernel.
    nan_policy: str
        How to handle NaN values. 'mask' will retain NaNs as NaNs, 'interp' will interpolate NaNs values.

    Returns
    -------
    np.ndarray
        The smoothed array.
    """
    if nan_policy == 'mask':
        truncate = 4.0
    elif nan_policy == 'interp':
        MAX_INTERP_DISTANCE = 100
        truncate = max(4.0, MAX_INTERP_DISTANCE/sigma)
    else:
        raise ValueError(f'Unknown nan_policy: {nan_policy}')

    v = x.copy()
    v[np.isnan(x)] = 0
    v_smooth = gaussian_filter(v, sigma=sigma, truncate=truncate)

    w = np.ones(x.shape)
    w[np.isnan(x)] = 0
    w_smooth = gaussian_filter(w, sigma=sigma, truncate=truncate)

    assert np.all((w_smooth > 0) | np.isnan(x)), 'Weights are zero at non-NaN positions.'
    w_smooth[w_smooth == 0] = 1

    out = v_smooth / w_smooth
    out[np.isnan(x)] = np.nan

    if nan_policy == 'mask':
        pass
    elif nan_policy == 'interp':
        out = interpolate_nans(out, sigma=sigma)
    else:
        raise ValueError(f'Unknown nan_policy: {nan_policy}')
    return out


def smoothing_sigma(k, scale=0.05, exponent=0.5, plateau=2):
    """
    Returns the standard deviation of the Gaussian kernel used for smoothing.

    The variance of the observed contact frequencies f are f(1-f), assuming a binomial distribution, yielding a
    coefficient of variation of sqrt(f(1-f))/f ~ 1/sqrt(f). (1 -f is approximately 1 for small f).

    We expect f ~ 1/k, where k is the genomic separation. Therefore, the coefficient of variation scales as
    1/sqrt(1/k) = sqrt(k).

    This motivates smoothing using a standard deviation proportional to sqrt(k).

    Since the contact frequencies do not decay significantly between k=1 and k=2, we set the standard deviation to
    a constant value for k < 2 (controlled by `plateau`).

    Parameters
    ----------
    k: int
        The distance between the two loci.
    scale: float
        The scaling factor for the standard deviation.
    exponent: float
        The exponent of the scaling factor.
    plateau: int
        The distance (in nucleosomes) where the scaling factor begins to increase.
    """
    return scale * np.maximum(1, k-plateau)**exponent


def coarsen_contact_map(contact_map, agg=np.nanmean, base=10, offset=0):
    """
    Smooth the contact map by aggregating the contact frequencies in blocks.

    Parameters
    ----------
    contact_map: np.ndarray
        The contact map to be coarsened.
    agg: callable
        The aggregation function.
    base: int
        Cells base**k from the diagonal are averaged over 2**k sized blocks.
    offset: int
        The offset from the diagonal where the first block starts.

    Returns
    -------
    np.ndarray
        The coarsened contact map.
    """
    assert base % 2 == 0

    n = contact_map.shape[0]
    result = contact_map.copy()
    max_k = int(np.ceil(np.log(n) / np.log(base)))

    for k in range(1, max_k):
        factor = 2 ** k  # How large a block to average over.
        start = offset + base ** k  # Offset from diagonal where first block starts
        end = offset + base ** (k+1)  # Offset from diagonal where the next blocks will start.
        for i in range(0, n - start, factor):
            for j in range(i+start, min(i+end, n), factor):
                result[i:i+factor, j:j+factor] = agg(contact_map[i:i+factor, j:j+factor])
    return np.triu(result, k=0) + np.triu(result, k=1).T


def smooth_diagonals(contact_map, scale=0.05, exponent=0.5, plateau=2, nan_policy='mask'):
    """
    Smooths the contact map using a Gaussian filter with a standard deviation that scales with the genomic separation.

    See `smoothing_sigma` for details on how the standard deviation is determined.

    Parameters
    ----------
    contact_map: np.ndarray
        The contact map to be smoothed.
    scale: float
        The scaling factor for the standard deviation.
    exponent: float
        The exponent of the scaling factor.
    plateau: int
        The distance (in nucleosomes) where the scaling factor begins to increase.
    nan_policy: str
        How to handle NaN values. 'mask' will retain NaNs as NaNs, 'interp' will interpolate NaNs values.

    Returns
    -------
    np.ndarray
        The smoothed contact map.
    """
    new_contact_map = np.zeros_like(contact_map)
    for k in range(contact_map.shape[0]):
        sigma = smoothing_sigma(k, scale=scale, exponent=exponent, plateau=plateau)
        smooth = nan_gaussian_filter(contact_map, sigma, nan_policy=nan_policy)
        new_diagonal = np.diag(smooth, k=k)
        for i in range(len(new_diagonal)):
            new_contact_map[i, i+k] = new_diagonal[i]
            new_contact_map[i+k, i] = new_diagonal[i]
    return new_contact_map


def interpolate_diagonals(contact_map, sigma=5):
    """
    Interpolate NaN values in the contact map by Gaussian smoothing along the diagonals.

    The motivation for smoothing along each diagonal is to preserve the P(s) curve as much as possible.

    Parameters
    ----------
    contact_map: np.ndarray
        The contact map to be smoothed.
    sigma: float or callable
        The standard deviation of the Gaussian filter.
    Returns
    -------
    np.ndarray
        The interpolated contact map.
    """
    new_contact_map = np.zeros_like(contact_map)
    ps = interpolate_nans(get_distance_average(contact_map, show_self=True), sigma=0)
    for k in range(contact_map.shape[0]):
        original_diagonal = np.diag(contact_map, k=k)
        if np.all(np.isnan(original_diagonal)):
            new_diagonal = np.zeros_like(original_diagonal) + ps[k]
        else:
            if isinstance(sigma, float) or isinstance(sigma, int):
                _sigma = sigma
            else:
                _sigma = sigma(k)
            new_diagonal = interpolate_nans(original_diagonal, sigma=_sigma)

        for i in range(len(new_diagonal)):
            new_contact_map[i, i+k] = new_diagonal[i]
            new_contact_map[i+k, i] = new_diagonal[i]
    return new_contact_map
