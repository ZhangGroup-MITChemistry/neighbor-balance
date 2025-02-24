"""
Test cases for experimental contact map balancing functions:

Balancing can be broken down into two steps:

1. Filtering: Removing bins that have low capture rates or low coverage.
2. Balancing: Adjusting the contact map to account for biases in the data.

Following filtering, it is expected that the contact map has no rows/columns that are all zeros. It is however,
reasonable to have rows/columns that are all NaNs, e.g. if a bin has very low coverage. This can be evaluated with the
`filtered_map_is_valid` function.

Key things to test related to filtering:
- min_nnz and mad_max filters are working as expected and respond to changes in the input parameters.
- These filters should be applied to both the rows and the columns of the contact map symmetrically.
- Capture rates should be combined in an additive way. I.e. if either bin has a non-zero capture rate, the capture rate
    for the joint bin should be non-zero.
- Having a zero joint capture rate for a bin should not cause all values in the contact map to be NaN.

Key things to test related to balancing:
- Entire diagonal is nan.
- Entire row is nan.
- A single entry is nan.
"""
import pytest
from neighbor_balance.ice import *


# Test get_marginal
example_maps = {
    'default': np.array([[0, 10, 5, 3],
                         [10, 0, 10, 5],
                         [5, 10, 0, 10],
                         [3, 5, 10, 0]]),
    'single_nan': np.array([[0, 10, np.nan, 3],
                            [10, 0, 10, 5],
                            [np.nan, 10, 0, 10],
                            [3, 5, 10, 0]]),
    'nan_row': np.array([[0, 10, np.nan, 3],
                         [10, 0, np.nan, 5],
                         [np.nan, np.nan, np.nan, np.nan],
                         [3, 5, np.nan, 0]]),
    'nan_diag': np.array([[0, 10, np.nan, 3],
                         [10, 0, 10, np.nan],
                         [np.nan, 10, 0, 10],
                         [3, np.nan, 10, 0]]),
    'nan_tip': np.array([[0, 10, 5, np.nan],
                          [10, 0, 10, 5],
                          [5, 10, 0, 10],
                          [np.nan, 5, 10, 0]]),
}


def test_get_marginal():
    contact_map = example_maps['default']
    marginal = get_marginal(contact_map)
    assert np.all(marginal == np.array([18, 25, 25, 18]))

    marginal = get_marginal(contact_map, k=2)
    assert np.all(marginal == np.array([8, 5, 5, 8]))


def test_get_marginal_single_nan():
    contact_map = example_maps['single_nan']
    marginal = get_marginal(contact_map)
    assert np.all(marginal == np.array([13, 25, 20, 18]))

    marginal = get_marginal(contact_map, k=2)
    assert np.all(marginal == np.array([3, 5, 0, 8]))


def test_get_marginal_single_nan_row():
    contact_map = example_maps['nan_row']
    marginal = get_marginal(contact_map)
    assert np.all(marginal == np.array([13, 15, 0, 8]))

    marginal = get_marginal(contact_map, k=2)
    assert np.all(marginal == np.array([3, 5, 0, 8]))


# Test interpolate_nans
def test_interpolate_nans_basic():
    x = np.array([0, 1, np.nan, 3, 4])
    x = interpolate_nans(x)
    assert np.all(~np.isnan(x))
    assert x[2] == 2


def test_interpolate_nans_all_nan():
    x = np.array([np.nan, np.nan, np.nan, np.nan, np.nan])
    with pytest.raises(ValueError):
        x = interpolate_nans(x)


def test_interpolate_nans_smoothing():
    x = np.array([2, 1, np.nan, np.nan, np.nan, np.nan, 3, 4])
    x = interpolate_nans(x, sigma=2)
    assert np.all(~np.isnan(x))
    assert np.allclose(x[2:6], np.linspace(1.5, 3.5, 6)[1:-1], atol=1e-1)


def test_interpolate_nans_no_smoothing():
    x = np.array([2, 1, np.nan, np.nan, np.nan, np.nan, 3, 4])
    x = interpolate_nans(x, sigma=0)
    assert np.all(~np.isnan(x))
    assert np.allclose(x[2:6], np.linspace(1, 3, 6)[1:-1], atol=1e-3)


# Test filter_bins
def test_filter_bins_min_nnz():
    contact_map = np.ones((100, 100))
    capture_rates = np.ones(contact_map.shape[0])
    contact_map = contact_map.astype(float)
    contact_map = filter_bins(contact_map, capture_rates, mad_max=0)
    assert filtered_map_is_valid(contact_map)
    assert np.all(~np.isnan(contact_map)), "All values should pass min_nnz"
    assert np.all(get_marginal(contact_map) == 99)

    contact_map = filter_bins(contact_map, capture_rates, min_nnz=200, mad_max=0)
    assert filtered_map_is_valid(contact_map)
    assert np.all(np.isnan(contact_map)), "All values should fail min_nnz and be NaN"
    assert np.all(get_marginal(contact_map) == 0)

    contact_map = np.ones((100, 100))
    contact_map[10] = 0
    contact_map[:, 10] = 0
    contact_map = contact_map.astype(float)
    contact_map = filter_bins(contact_map, capture_rates, min_nnz=90, mad_max=0)

    assert filtered_map_is_valid(contact_map)
    assert np.all(np.isnan(contact_map[10])), "Row 10 should fail min_nnz and be NaN"
    assert np.all(np.isnan(contact_map[:, 10]))
    contact_map[10] = 0
    contact_map[:, 10] = 0
    assert np.all(~np.isnan(contact_map)), "All other values should pass min_nnz"


def test_filter_bins_mad_max():
    # TODO this is a bad test because the mad is zero, so we can't test that the threshold is working.
    contact_map = np.ones((100, 100))
    contact_map[10] = 0.1
    contact_map[:, 10] = 0.1

    capture_rates = np.ones(contact_map.shape[0])
    contact_map = contact_map.astype(float)
    contact_map = filter_bins(contact_map, capture_rates)
    assert filtered_map_is_valid(contact_map)
    assert np.all(np.isnan(contact_map[10])), "Row 10 should fail mad_max and be NaN"
    assert np.all(np.isnan(contact_map[:, 10]))
    contact_map[10] = 0
    contact_map[:, 10] = 0
    assert np.all(~np.isnan(contact_map)), "All other values should pass mad_max"


def test_filter_bins_mad_max_off():
    contact_map = np.ones((100, 100))
    contact_map[10] = 0.1
    contact_map[:, 10] = 0.1

    capture_rates = np.ones(contact_map.shape[0])
    contact_map = contact_map.astype(float)
    contact_map = filter_bins(contact_map, capture_rates, mad_max=0)
    assert filtered_map_is_valid(contact_map)
    assert np.all(~np.isnan(contact_map)), "All other values should pass mad_max"


def test_filter_bins_capture_rates_one_zero():
    contact_map = np.ones((100, 100))
    capture_rates = np.ones(contact_map.shape[0])
    capture_rates[10] = 0
    contact_map = contact_map.astype(float)
    contact_map = filter_bins(contact_map, capture_rates)
    assert filtered_map_is_valid(contact_map)
    contact_map[10, 10] = 0  # Only the self-interaction should be zero, but this doesn't matter.
    assert np.all(~np.isnan(contact_map)), "Having one zero capture rate should not cause all values to be NaN"


def test_filter_bins_capture_rates_two_zero():
    contact_map = np.ones((100, 100))
    capture_rates = np.ones(contact_map.shape[0])
    capture_rates[10] = 0
    capture_rates[20] = 0
    contact_map = contact_map.astype(float)
    contact_map = filter_bins(contact_map, capture_rates, mad_max=0)
    assert filtered_map_is_valid(contact_map)
    contact_map[10, 10] = 0  # Only the self-interaction should be zero, but this doesn't matter.
    contact_map[20, 20] = 0
    assert np.isnan(contact_map[10, 20]), "Having two zero capture rate should cause one NaN"
    assert np.isnan(contact_map[20, 10])
    contact_map[10, 20] = 0
    contact_map[20, 10] = 0
    assert np.all(~np.isnan(contact_map)), "Having two zero capture rate should not cause all values to be NaN"


def test_filter_bins_large_chunk_of_zero_capture_rates():
    contact_map = np.ones((100, 100))
    capture_rates = np.ones(contact_map.shape[0])
    capture_rates[:50] = 0
    contact_map = contact_map.astype(float)
    contact_map = filter_bins(contact_map, capture_rates, min_nnz=40, mad_max=0)
    assert filtered_map_is_valid(contact_map)

    assert np.all(np.isnan(contact_map[:50, :50])), "All entries with low joint capture rate should be NaN"
    assert np.all(~np.isnan(contact_map[50:, 50:])), "All entries with high joint capture rate should be non-NaN"
    assert np.all(~np.isnan(contact_map[:50, 50:])), "All entries with low marginal capture rate should be non-NaN"
    assert np.all(~np.isnan(contact_map[50:, :50])), "All entries with low marginal capture rate should be non-NaN"


# We get different results if filters are applied in different orders.
# The order I've chosen makes this test obsolete.
# def test_filter_bins_large_chunk_of_zero_capture_rates_puts_below_nnz():
#     contact_map = np.ones((100, 100))
#     capture_rates = np.ones(contact_map.shape[0])
#     capture_rates[:50] = 0
#     contact_map = contact_map.astype(float)
#     contact_map = filter_bins(contact_map, capture_rates, min_nnz=60, mad_max=0)
#     assert filtered_map_is_valid(contact_map)
#     assert np.all(np.isnan(contact_map)), "All entries with low joint capture rate should be NaN"


# Test ice_balance

def test_ice_balance():
    contact_map = example_maps['default'].astype(float)
    assert filtered_map_is_valid(contact_map)
    contact_map, _ = ice_balance(contact_map)
    assert np.all(~np.isnan(contact_map))


def test_ice_balance_single_nan():
    contact_map = example_maps['single_nan'].astype(float)
    assert filtered_map_is_valid(contact_map)
    contact_map, _ = ice_balance(contact_map)
    assert np.isnan(contact_map[0, 2]), "NaN values should remain NaN"
    contact_map[0, 2] = 0
    contact_map[2, 0] = 0
    assert np.all(~np.isnan(contact_map)), 'All values should be non-NaN'


def test_ice_balance_row_nan():
    contact_map = example_maps['nan_row'].astype(float)
    assert filtered_map_is_valid(contact_map)
    contact_map, _ = ice_balance(contact_map)
    assert np.all(np.isnan(contact_map[2])), "NaN values should remain NaN"
    assert np.all(np.isnan(contact_map[:, 2]))
    contact_map[2] = 0
    contact_map[:, 2] = 0
    assert np.all(~np.isnan(contact_map)), 'All values should be non-NaN'


def test_ice_balance_interpolate():
    contact_map = example_maps['single_nan'].astype(float)
    contact_map, _ = ice_balance(contact_map, interpolate=True)
    assert np.all(~np.isnan(contact_map))


def test_ice_balance_interpolate_nan_diag():
    contact_map = example_maps['nan_diag'].astype(float)
    contact_map, _ = ice_balance(contact_map, interpolate=True)
    assert np.all(~np.isnan(contact_map))


def test_ice_balance_interpolate_nan_tip():
    contact_map = example_maps['nan_tip'].astype(float)
    contact_map, _ = ice_balance(contact_map, interpolate=True)
    assert np.all(~np.isnan(contact_map))


def test_ice_balance_zeros():
    contact_map = np.array([[0, 0, 0, 0],
                            [0, 0, 10, 5],
                            [0, 10, 0, 10],
                            [0, 5, 10, 0]])
    with pytest.raises(ValueError):
        contact_map, _ = ice_balance(contact_map)
