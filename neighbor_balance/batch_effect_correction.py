import os
import numpy as np
import click

from neighbor_balance.plotting import ContactMap
from neighbor_balance.neighbor import normalize_contact_map_neighbor


def get_ps(contact_map, average=np.mean):
    n = contact_map.shape[0]
    ps = np.zeros(n)
    for i in range(n):
        ps[i] = average(np.diagonal(contact_map, offset=i))
    return ps


def get_average_ps(contact_maps, high=None):
    ps = []
    for contact_map in contact_maps:
        ps.append(get_ps(contact_map)[:high])
    return np.mean(ps, axis=0)


def get_scaling_factors(query_contact_maps, reference_contact_maps, reference_separation=5):
    query_ps = get_average_ps(query_contact_maps, high=reference_separation + 1)
    reference_ps = get_average_ps(reference_contact_maps, high=reference_separation + 1)
    global_factor = np.mean(reference_ps[reference_separation]) / np.mean(query_ps[reference_separation])
    scaled_ps = query_ps * global_factor

    factors = []
    for separation in range(reference_separation):
        factors += [reference_ps[separation] / scaled_ps[separation]]
    return global_factor, factors


def remove_batch_effects(query_contact_map, global_factor, factors):
    scaled_contact_map = query_contact_map * global_factor
    seps = np.abs(
        np.arange(scaled_contact_map.shape[0])[:, None]
        - np.arange(scaled_contact_map.shape[1])[None, :]
    )
    for separation, factor in enumerate(factors):
        scaled_contact_map[seps == separation] *= factor
    return scaled_contact_map


def compute_factors(contact_maps: dict, reference: str, query: str, regions: list, reference_separation: int = 5):
    """Hi-C: midG1 reference, anatelo query — per-region and merged."""
    factors_dict = {}
    for region in regions:
        reference_cmap = contact_maps[region][reference].contact_map.copy()
        query_cmap = contact_maps[region][query].contact_map.copy()
        global_factor, factors = get_scaling_factors(
            [query_cmap],
            [reference_cmap],
            reference_separation=reference_separation,
        )
        factors_dict[region] = {
            "global_factor": global_factor,
            "factors": factors,
        }
    all_reference_cmaps = [contact_maps[region][reference].contact_map for region in regions]
    all_query_cmaps = [contact_maps[region][query].contact_map for region in regions]
    global_factor, factors = get_scaling_factors(
        all_query_cmaps,
        all_reference_cmaps,
        reference_separation=reference_separation,
    )
    factors_dict["merged"] = {
        "global_factor": global_factor,
        "factors": factors,
    }
    return factors_dict


def print_factors(factors_dict):
    for region, factors_info in factors_dict.items():
        global_factor = factors_info["global_factor"]
        factors = factors_info["factors"]
        print(f"region: {region}")
        print(f"global_factor: {global_factor:.3f}")
        print(f"factors: {','.join(f'{factor:.3f}' for factor in factors)}")


def load_contact_maps(conditions, regions, root_dir, pattern="{root_dir}/{condition}/{condition}_output_{region}.npz"):
    contact_maps = {}
    for region in regions:
        contact_maps[region] = {}
        for condition in conditions:
            print(region, condition)
            contact_map_file = pattern.format(root_dir=root_dir, condition=condition, region=region)
            cmap = ContactMap.from_npz(contact_map_file)
            cmap.contact_map = normalize_contact_map_neighbor(cmap.contact_map)
            contact_maps[region][condition] = cmap
    return contact_maps


@click.command()
@click.option('--reference', default='midG1', help='Reference condition for scaling factors.')
@click.option('--root', default='/orcd/data/binz/001/belong/mouse_cell_cycle_PIPELINE', help='Root directory for contact map files.')
@click.option('--reference-separation', default=5, help='Separation distance for reference scaling factors.')
@click.option('--conditions', default='async,prometa,anatelo,earlyG1,midG1,lateG1,smc2-aid-0hr,smc2-aid-4hr', help='Comma-separated list of conditions.')
@click.option('--regions', default='1,2,3,4,5', help='Comma-separated list of regions.')
@click.option('--pattern', default='{root_dir}/{condition}/{condition}_output_{region}.npz', help='Pattern for contact map file paths.')
def main(reference, root, reference_separation, conditions, regions, pattern):
    conditions = conditions.split(',')
    regions = [r for r in regions.split(',')]

    contact_maps = load_contact_maps(
        conditions=conditions,
        regions=regions,
        root_dir=root,
        pattern=pattern
    )

    for query in conditions:
        print(f"##### query={query} reference={reference}")
        factors_dict = compute_factors(contact_maps, reference, query, regions=regions, reference_separation=reference_separation)
        print_factors(factors_dict)

if __name__ == "__main__":
    main()
