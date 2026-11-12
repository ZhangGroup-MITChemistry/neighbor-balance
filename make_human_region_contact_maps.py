import os
import numpy as np
import matplotlib.pyplot as plt

from neighbor_balance.plotting import parse_region, ContactMap
from neighbor_balance.smoothing import coarsen_contact_map
from neighbor_balance.neighbor import normalize_contact_map_neighbor, normalize_contact_map_average


condition = 'GM12878'
output_directory = '.'

# Available on Zenodo, update input directory to the location where the data is downloaded.
input_directory = f"/home/joepaggi/binz/belong/for_revisions"
input_cool = f"{input_directory}/{condition}_minus_inward.mcool" 

# See Supplementary Table 3 from Hong 2025, "Genome structure mapping with high-resolution 3D genomics and deep learning".
capture_probes_bed = "/orcd/data/binz/001/joepaggi/rcmc_from_clarice/supp_table3_regions_covered_by_probes.bed"

regions = {
  '1': "chr6:25,112,000-28,621,800",
  '3': "chr5:157,036,200-160,147,400",
  '4': "chr1:207,626,400-210,338,200",
  '5': "chr4:61,369,200-64,434,200",
  '6': "chr8:124,988,400-129,795,600",
  '7': "chr6:29,678,600-32,256,200",
  '8': "chrX:47,081,600-49,442,000",
  '9': "chr1:237,287,000-240,531,000",
  '10': "chr7:10,480,200-13,506,200",
  '11': "chr8:62,919,000-66,566,200",
  '12': "chr4:181,290,000-184,014,400",
  '13': "chr3:118,624,600-121,534,000",
  '14': "chr9:106,603,400-109,909,200",
}

print('Submit these separately for each region with your job scheduler...')
for name, region in regions.items():
    out_file = f"{output_directory}/{condition}_{name}.npz"
    cmd = f"neighbor-balance region-balance --capture-probes-path {capture_probes_bed} {out_file} {input_cool} {region}"
    print(cmd)


# Plot the contact density for an example region.
name = '1'
out_file = f"{output_directory}/{condition}_{name}.npz"
reference_out_file = f'/orcd/data/binz/001/belong/for_revisions/{condition}/{condition}/{condition}_output_{name}.npz' # Used for figures in paper.

if not os.path.exists(out_file):
    print(f"# {out_file} does not exist.")
    exit()


def load_and_process_region(out_file):
    loaded = np.load(out_file, allow_pickle=True)
    contact_map = loaded['array']  # Unnormalized ICE balanced contact map as nxn numpy array.
    metadata = loaded['metadata'].item()
    chrom, start, end = parse_region(metadata['region'])

    ice_contact_map = normalize_contact_map_average(contact_map.copy())
    ice_contact_map = coarsen_contact_map(ice_contact_map)
    ice_contact_map = ContactMap(ice_contact_map, chrom, start, end, 200)

    neighbor_contact_map = normalize_contact_map_neighbor(contact_map.copy())
    neighbor_contact_map = coarsen_contact_map(neighbor_contact_map)
    neighbor_contact_map = ContactMap(neighbor_contact_map, chrom, start, end, 200)
    return neighbor_contact_map, ice_contact_map


neighbor_contact_map, ice_contact_map = load_and_process_region(out_file)
f, ax = plt.subplots(figsize=(7, 2))
ax.plot(neighbor_contact_map.get_marginal())

# Show the marginal of the reference file as a sanity check.
if os.path.exists(reference_out_file):
    reference_neighbor_contact_map, reference_ice_contact_map = load_and_process_region(reference_out_file)
    ax.plot(reference_neighbor_contact_map.get_marginal(), label='reference', ls='--')
  
plt.savefig(f"{output_directory}/{condition}_{name}_marginal.png")
plt.close()
