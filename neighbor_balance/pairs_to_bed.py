"""
Example usage:

chromsizes=~/ucsc/mm39.chrom.sizes
sbatch -J mnase --cpus-per-task 10 \
    --wrap="zcat ../WT_BR1/all.nodups.pairs.gz |
    python ~/na_genes/na_genes/microc/pairs_to_bed.py $chromsizes --protected-over-2 73 | \
    sort -k1,1 -k2,2n | \
    bedtools genomecov -i - -g $chromsizes -5 -bg > all.nodups.shift73.bg"
sbatch -J mnase --dependency=SINGLETON \
    --wrap="bedGraphToBigWig all.nodups.shift73.bg $chromsizes all.nodups.shift73.bw"

Can add `grep +` before sort to filter for reads that are on the plus strand
or `grep -` for reads that are on the minus strand.
"""

import click
import sys
import pandas as pd


@click.command()
@click.argument('chromsizes')
@click.option('--i', help='Input file in pairs format, stdin if not specified.', default=None)
@click.option('--o', help='Output file in bed format, stdout if not specified.', default=None)
@click.option('--protected-over-2', help='Shift each position 3\' by this amount', default=0)
def main(chromsizes, i, o, protected_over_2):
    """
    Convert pairs format file to bed format file with one entry for each read end.
    """
    chromsizes = pd.read_csv(chromsizes, sep='\t', header=None, names=['chrom', 'size'])
    chromsizes = {k: v for k, v in zip(chromsizes['chrom'], chromsizes['size'])}

    if i is None:
        i = sys.stdin
    else:
        i = open(i, 'r')

    if o is None:
        o = sys.stdout
    else:
        o = open(o, 'w')

    for line in i:
        if line[0] == '#':
            continue

        name, chrom1, pos1, chrom2, pos2, strand1, strand2, read_type = line.strip().split('\t')[:8]
        if read_type != 'UU':
            continue

        if strand1 == '+':
            pos1 = int(pos1) + protected_over_2
        else:
            pos1 = int(pos1) - protected_over_2
        if strand2 == '+':
            pos2 = int(pos2) + protected_over_2
        else:
            pos2 = int(pos2) - protected_over_2

        if pos1 < 0 or pos2 < 0 or pos1 >= chromsizes[chrom1] or pos2 >= chromsizes[chrom2]:
            continue

        o.write(f'{chrom1}\t{pos1}\t{pos1+1}\t{strand1}\n')
        o.write(f'{chrom2}\t{pos2}\t{pos2+1}\t{strand2}\n')

    i.close()
    o.close()


if __name__ == '__main__':
    main()
