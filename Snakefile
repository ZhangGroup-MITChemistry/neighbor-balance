condition_done_files = [
    f"{condition}/{condition}_minus_inward_neighbor.done"
    for condition in config["samples"].keys()
]

replicate_done_files = [
    f"{condition}/{replicate}/{condition}_{replicate}_minus_inward_neighbor.done"
    for condition in config["samples"].keys()
    for replicate in config["samples"][condition].keys()
]

base_ps_files = [
    f"{condition}/{replicate}/{run}/{run}.pairs.gz.base_ps.csv"
    for condition in config["samples"].keys()
    for replicate in config["samples"][condition].keys()
    for run in config["samples"][condition][replicate]
]

rule all:
    input:
        condition_done_files,
        replicate_done_files,
        base_ps_files

################################################################################################################
# Per-run rules (different files for the same experiment)

rule download_fastq:
    output:
        fastq1='{condition}/{replicate}/{run}/{run}_1.fastq',
        fastq2='{condition}/{replicate}/{run}/{run}_2.fastq'
    resources:
        mem_mb=32000
    shell:
        """
        mkdir -p {wildcards.condition}/{wildcards.replicate}/{wildcards.run}
        fasterq-dump {wildcards.run} -O {wildcards.condition}/{wildcards.replicate}/{wildcards.run}
        """

rule make_pairs:
    input:
        fastq1='{run}_1.fastq',
        fastq2='{run}_2.fastq',
    output:
        pairs="{run}.pairs.gz"
    resources:
        mem_mb=32000
    threads: 48
    shell:
        """
        tmpfile=$(mktemp {output.pairs}.tmp.XXXXXX.gz)
        bowtie2 -x {config[index]} --threads {threads} -1 {input.fastq1} -2 {input.fastq2} --reorder --local --very-sensitive-local --maxins 1 --minins 1000000 \
         | pairtools parse --add-columns mapq --walks-policy mask -c {config[chromsizes]} --assembly {config[genome]} --min-mapq 2 --drop-sam \
         | pairtools sort -o $tmpfile
        mv $tmpfile {output.pairs}
        """

rule base_ps:
    input:
        pairs="{run}.pairs.gz"
    output:
        base_ps="{run}.pairs.gz.base_ps.csv"
    resources:
        mem_mb=32000
    shell:
        """
        neighbor-balance base-ps {input.pairs} --output {output.base_ps}
        """

################################################################################################################
# Per-experiment rules
rule merge_and_process_pairs:
    input:
        pairs=lambda wildcards: expand(
            "{condition}/{replicate}/{run}/{run}.pairs.gz",
            condition=wildcards.condition,
            replicate=wildcards.replicate,
            run=config["samples"][wildcards.condition][wildcards.replicate]
        ),
        chromsizes=config["chromsizes"]
    output:
        merged_pairs="{condition}/{replicate}/{condition}_{replicate}.nodups.select.shifted.pairs.gz",
        stats="{condition}/{replicate}/{condition}_{replicate}.dedup.stats"
    threads: 8
    resources:
        mem_mb=32000
    params:
        capture_opt = f"--regions {config['regions']}" if "regions" in config else ""
    shell:
        """
        tmp_pairs=$(mktemp {output.merged_pairs}.tmp.XXXXXX.gz)
        tmp_stats=$(mktemp {output.stats}.tmp.XXXXXX)
        pairtools merge {input.pairs} \
        | pairtools dedup --max-mismatch 1 --mark-dups --output-stats $tmp_stats \
        | neighbor-balance shift-pairs --protected-over-2 65 \
        | neighbor-balance filter-pairs {input.chromsizes} {params.capture_opt}\
        | pairtools sort -o $tmp_pairs
        mv $tmp_pairs {output.merged_pairs}
        mv $tmp_stats {output.stats}
        """

rule generate_cool:
    input:
        merged_pairs="{condition}/{replicate}/{condition}_{replicate}.nodups.select.shifted.pairs.gz",
        chromsizes=config["chromsizes"]
    output:
        cool="{condition}/{replicate}/{condition}_{replicate}.cool"
    threads: 6
    resources:
        mem_mb=32000
    shell:
        """
        zcat {input.merged_pairs} \
            | cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 {input.chromsizes}:200 - {output.cool}
        """

rule generate_cool_inward:
    input:
        merged_pairs="{condition}/{replicate}/{condition}_{replicate}.nodups.select.shifted.pairs.gz",
        chromsizes=config["chromsizes"]
    output:
        cool_inward="{condition}/{replicate}/{condition}_{replicate}_inward.cool"
    threads: 6
    resources:
        mem_mb=32000
    shell:
        """
        zcat {input.merged_pairs} \
            | neighbor-balance filter-pairs {input.chromsizes} --direction inward \
            | cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 {input.chromsizes}:200 - {output.cool_inward}
        """

rule remove_inward:
    input:
        cool="{condition}/{replicate}/{condition}_{replicate}.cool",
        cool_inward="{condition}/{replicate}/{condition}_{replicate}_inward.cool"
    output:
        cool_minus_inward="{condition}/{replicate}/{condition}_{replicate}_minus_inward.cool"
    resources:
        mem_mb=32000
    shell:
        """
        neighbor-balance remove-inward {input.cool} {input.cool_inward} {output.cool_minus_inward}
        """

################################################################################################################
# Per-condition rules
rule merge_cools:
    input:
        cools=lambda wildcards: expand(
            "{condition}/{replicate}/{condition}_{replicate}_minus_inward.cool",
            condition=wildcards.condition,
            replicate=config["samples"][wildcards.condition].keys()
        ),
    output:
        merged="{condition}/{condition}_minus_inward.cool"
    resources:
        mem_mb=32000
    shell:
        """
        cooler merge {output.merged} {input.cools}
        """

rule zoomify:
    input:
        "{name}_minus_inward.cool"
    output:
        "{name}_minus_inward.mcool"
    threads: 8
    resources:
        mem_mb=32000
    shell:
        """
        cooler zoomify --balance --balance-args '--ignore-diags 0' \
            --resolutions 200,400,800,1600,3200,6400,12800,25600,51200,102400 \
            --out {output} {input}
        """

rule neighbor_balance:
    input:
        "{name}_minus_inward.mcool"
    output:
        "{name}_minus_inward_neighbor.done"
    resources:
        mem_mb=32000
    shell:
        """
        neighbor-balance neighbor-balance-cooler {input}
        touch {output}
        """
