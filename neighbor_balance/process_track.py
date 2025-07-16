import click
import sys
import gzip
import os

def _select_wig_track(fp, track):
    """
    Select a specific track from a WIG file and write it to stdout.

    Parameters
    ----------
    fp : file-like object
        File-like object to read from.
    track : str
        Name of the track to select.
    """
    detected_track = False
    for line in fp:
        if detected_track:
            if line.startswith("track"):
                break
            print(line.strip())
        elif line.startswith("track"):
            detected_track = (track is None) or (track in line)


def _bedgraph_merge(fp):
    """
    Merge overlapping intervals in a bedgraph file to facilitate conversion to bigwig.
    """
    prev_chr = None
    prev_chr_e = None
    for line in fp:
        fields = line.strip().split('\t')

        if len(fields) < 3:
            raise ValueError("Invalid bedgraph line: %s" % line)

        chr_name = fields[0]
        chr_start = int(fields[1])
        chr_end = int(fields[2])

        if prev_chr is None or prev_chr != chr_name or prev_chr_e <= chr_start:
            print(line.strip())

            prev_chr = chr_name
            prev_chr_e = chr_end


@click.group()
def cli():
    pass


@cli.command()
@click.argument('filename', type=click.Path(exists=True))
@click.option('--track', default=None, help='Name of the track to select. If not provided, the first track is selected.')
def select_wig_track(filename, track):
    if filename == "-":
        _select_wig_track(sys.stdin, track)
    elif filename.endswith(".gz"):
        with gzip.open(filename, 'rt') as fp:
            _select_wig_track(fp, track)
    else:
        with open(filename, 'r') as fp:
            _select_wig_track(fp, track)


@cli.command()
@click.argument('filename', type=click.Path(exists=True))
def bedgraph_merge(filename):
    if filename == "-":
        _bedgraph_merge(sys.stdin)
    else:
        with open(filename, 'r') as fp:
            _bedgraph_merge(fp)


@cli.command()
@click.argument('name')
@click.argument('ftp_url')
@click.option('--input-type', default='bw')
@click.option('--input-genome', default='mm10')
@click.option('--missing-chr', is_flag=True)
@click.option('--no-track', is_flag=True)
@click.option('--chain-dir', default='~/ucsc', help='Directory containing chain files for liftover.')
def convert(name, ftp_url, input_type, input_genome, missing_chr, no_track, chain_dir):
    assert input_type in ['bw', 'wig.gz', 'bg.gz']

    # Download the file. Compute nodes don't have internet, so we need to download it on the login node.
    print(f'wget -O temp.{name}.{input_type} {ftp_url}')

    cmds = []
    # Convert to bedgraph with no track lines.
    if input_type == 'wig.gz':
        if no_track:
            cmds += [f'zcat temp.{name}.wig.gz > temp.{name}.wig']
        else:
            cmds += [f'python {__file__} select-wig-track temp.{name}.wig.gz > temp.{name}.wig']
        cmds += [
            f'wigToBigWig temp.{name}.wig  {chain_dir}/{input_genome}.chrom.sizes temp.{name}.bw -clip',
            f'bigWigToBedGraph temp.{name}.bw temp.{name}.bg'
        ]
    elif input_type == 'bg.gz':
        cmds += [f'python {__file__} select-wig-track temp.{name}.bg.gz > temp.{name}.bg']
    else:
        cmds += [f'bigWigToBedGraph temp.{name}.bw temp.{name}.bg']

    if missing_chr:
        cmds += [f"sed -i -e 's/^/chr/' temp.{name}.bg"]

    # Lift over to mm39.
    if input_genome == 'mm8':
        cmds += [
            f'liftOver temp.{name}.bg {chain_dir}/mm8ToMm10.over.chain.gz temp.{name}.mm10.bg temp.{name}.mm10.unmap.bg',
            f'liftOver temp.{name}.bg {chain_dir}/mm10ToMm39.over.chain.gz temp.{name}.mm39.bg temp.{name}.mm39.unmap.bg',
        ]
    elif input_genome in ['mm9', 'mm10']:
        cmds += [
            f'liftOver temp.{name}.bg {chain_dir}/{input_genome}ToMm39.over.chain.gz temp.{name}.mm39.bg temp.{name}.mm39.unmap.bg',
        ]
    else:
        raise ValueError(f'Invalid input genome: {input_genome}')

    # Convert to bigwig.
    cmds += [
        f'sort -k1,1 -k2,2n temp.{name}.mm39.bg > temp.{name}.mm39.sorted.bg',
        f'python {__file__} bedgraph-merge temp.{name}.mm39.sorted.bg >  temp.{name}.mm39.sorted.merged.bg',
        f'bedGraphToBigWig temp.{name}.mm39.sorted.merged.bg {chain_dir}/mm39.chrom.sizes {name}.mm39.bw',
    ]

    cmds = ';'.join(cmds)
    assert '"' not in cmds, "Double quotes are not allowed in the command string."
    print(f'sbatch -J {name} -o {name}.log --wrap="{cmds}"')


@cli.command()
@click.argument('fname')
def mm39tomm10(fname):
    assert fname.endswith('.mm39.bw')
    name = fname[:-8]

    cmds = []

    # bw to bg
    cmds += [f'bigWigToBedGraph {name}.mm39.bw temp.{name}.mm39.bg']

    # liftover
    cmds += [f'liftOver temp.{name}.mm39.bg {chain_dir}/mm39ToMm10.over.chain.gz temp.{name}.mm10.bg temp.{name}.mm10.unmap.bg']

    # bg to bw
    cmds += [
        f'sort -k1,1 -k2,2n temp.{name}.mm10.bg > temp.{name}.mm10.sorted.bg',
        f'python {__file__} bedgraph-merge temp.{name}.mm10.sorted.bg >  temp.{name}.mm10.sorted.merged.bg',
        f'bedGraphToBigWig temp.{name}.mm10.sorted.merged.bg {chain_dir}/mm10.chrom.sizes {name}.mm10.bw',
    ]

    cmds = ';'.join(cmds)
    assert '"' not in cmds, "Double quotes are not allowed in the command string."
    print(f'sbatch -J {name} -o {name}.mm39tomm10.log --wrap="{cmds}"')


if __name__ == "__main__":
    cli()
