"""
Convert a pairs file to BED format using median-assigned positions.

This assigns each read end to the median of its mapped 5' and 3' positions
when available (pos51/pos31 for read1 and pos52/pos32 for read2). If those
columns are missing or set to '.', falls back to the existing pos1/pos2.

Example usage:

chromsizes=~/ucsc/mm39.chrom.sizes
zcat ../WT_BR1/all.nodups.pairs.gz \
| pairs_to_bed_median $chromsizes \
| sort -k1,1 -k2,2n \
| bedtools genomecov -i - -g $chromsizes -5 -bg \
> all.nodups.median.bg

bedGraphToBigWig all.nodups.median.bg $chromsizes all.nodups.median.bw
"""

import sys
import click
import pandas as pd

from .pairs import PairsFile, shift_line_median


@click.command()
@click.argument('chromsizes')
@click.option('--i', help='Input file in pairs format, stdin if not specified.', default=None)
@click.option('--o', help='Output file in BED format, stdout if not specified. Ignored if --length-threshold is set and split outputs are provided.', default=None)
@click.option('--pair-type', default='UU', help="Only output pairs with this pair_type (set empty to disable).")
@click.option('--length-threshold', type=int, default=None, help='Split outputs by abs(pos5-pos3) per end. Values < threshold go to BELOW, >= threshold to ABOVE.')
@click.option('--o-below', default=None, help='Output BED for ends with length < threshold (required if --length-threshold is set).')
@click.option('--o-above', default=None, help='Output BED for ends with length >= threshold (required if --length-threshold is set).')
def main(chromsizes, i, o, pair_type, length_threshold, o_below, o_above):
    """
    Convert pairs format file to BED format using median-assigned positions.
    Writes one 1bp BED entry per read end.
    """
    chromsizes_df = pd.read_csv(chromsizes, sep='\t', header=None, names=['chrom', 'size'])
    chromsizes_dict = {k: v for k, v in zip(chromsizes_df['chrom'], chromsizes_df['size'])}

    # Input
    if i is None:
        fp_in = sys.stdin
        close_in = False
    else:
        fp_in = open(i, 'r')
        close_in = True

    # Outputs
    if length_threshold is None:
        # Single-stream output (legacy behavior)
        if o is None:
            fp_out = sys.stdout
            close_out = False
        else:
            fp_out = open(o, 'w')
            close_out = True
        fp_out_below = None
        fp_out_above = None
    else:
        # Split outputs by threshold; both outputs required
        if o_below is None or o_above is None:
            raise click.UsageError('When --length-threshold is set, both --o-below and --o-above must be provided.')
        fp_out_below = open(o_below, 'w')
        fp_out_above = open(o_above, 'w')
        fp_out = None
        close_out = False

    pf = PairsFile(fp_in)

    # Helper to parse integer position or None
    def _parse(value):
        return None if value == '.' else int(value)

    for line in pf:
        # Optional pair_type filter (e.g., keep only uniquely mapping reads 'UU')
        if pair_type and line.get('pair_type', None) != pair_type:
            continue

        # Apply median assignment
        line = shift_line_median(line)

        # Extract fields
        chrom1 = line['chrom1']
        chrom2 = line['chrom2']
        strand1 = line['strand1']
        strand2 = line['strand2']
        pos1 = line.get('pos1', '.')
        pos2 = line.get('pos2', '.')

        # Skip if positions are missing
        if pos1 == '.' or pos2 == '.':
            continue

        pos1 = int(pos1)
        pos2 = int(pos2)

        # Bounds check
        if chrom1 not in chromsizes_dict or chrom2 not in chromsizes_dict:
            continue
        if pos1 < 0 or pos1 >= chromsizes_dict[chrom1]:
            continue
        if pos2 < 0 or pos2 >= chromsizes_dict[chrom2]:
            continue

        if length_threshold is None:
            # Single-stream output
            fp_out.write(f"{chrom1}\t{pos1}\t{pos1+1}\t{strand1}\n")
            fp_out.write(f"{chrom2}\t{pos2}\t{pos2+1}\t{strand2}\n")
        else:
            # Compute end-specific alignment lengths when possible
            len1 = None
            if 'pos51' in line and 'pos31' in line:
                p51 = _parse(line['pos51'])
                p31 = _parse(line['pos31'])
                if p51 is not None and p31 is not None:
                    len1 = abs(p51 - p31)

            len2 = None
            if 'pos52' in line and 'pos32' in line:
                p52 = _parse(line['pos52'])
                p32 = _parse(line['pos32'])
                if p52 is not None and p32 is not None:
                    len2 = abs(p52 - p32)

            # Route each end to appropriate file if its length is known
            if len1 is not None:
                if len1 < length_threshold:
                    fp_out_below.write(f"{chrom1}\t{pos1}\t{pos1+1}\t{strand1}\n")
                else:
                    fp_out_above.write(f"{chrom1}\t{pos1}\t{pos1+1}\t{strand1}\n")

            if len2 is not None:
                if len2 < length_threshold:
                    fp_out_below.write(f"{chrom2}\t{pos2}\t{pos2+1}\t{strand2}\n")
                else:
                    fp_out_above.write(f"{chrom2}\t{pos2}\t{pos2+1}\t{strand2}\n")

    if close_in:
        fp_in.close()
    if close_out:
        fp_out.close()
    if length_threshold is not None:
        fp_out_below.close()
        fp_out_above.close()


if __name__ == '__main__':
    main()




