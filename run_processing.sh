echo "This is not meant to be run as a script. Please run commands individually."
exit

####################################################################################################
# Define variables and
conda activate neighbor
genome=mm39
index=~/ucsc/mm39/mm39
chromsizes=~/ucsc/mm39.chrom.sizes
captureregions=~/rcmc/captureregions_mm39.bed
root=~/rcmc/WT_BR1
name=WT_BR1


###################################################################################################
# Download raw fastq files data from the SRA and capture pairs file from github.

cd $root || exit
wget https://raw.githubusercontent.com/ahansenlab/RCMC_analysis_code/refs/heads/main/captureprobes_mm39.bed
prefetch -p --max-size 100GB SRX15946369 # WT BR1


####################################################################################################
# Unpack SRA archives.
cd "$root" || exit
for sra in SRR*; do
    cd $root/$sra || exit
    sbatch -J $sra -p $partition --dependency=SINGLETON --wrap="fasterq-dump -p ${sra}.sra"
done


####################################################################################################
# Align and convert to pairs.
cd "$root" || exit
for sra in SRR*; do
    cd "$root/$sra" || exit
    cmd="bowtie2 -x $index --threads 48 -1 ${sra}_1.fastq -2 ${sra}_2.fastq --reorder --local --very-sensitive-local --maxins 1 --minins 1000000 \
         | pairtools parse --add-columns mapq --walks-policy mask -c $chromsizes --assembly $genome --min-mapq 2 --drop-sam \
         | pairtools sort -o ${sra}.pairs.gz"
    echo "$cmd"
    sbatch --cpus-per-task 48 -p $partition -J $sra --dependency=SINGLETON --wrap="$cmd"
done


####################################################################################################
# Merge, dedup, shift, and filter.
# For genome-wide datasets, remove the --regions argument.
sbatch  --mem=32G--cpus-per-task 16 -p $partition -t 96:00:00 -J $name --dependency=SINGLETON -p $partition \
  --wrap="pairtools merge SRR*/SRR*.pairs.gz \
  | pairtools dedup --max-mismatch 1 --mark-dups --output-stats ${name}.dedup.stats \
  | neighbor-balance shift-pairs --protected-over-2 65 \
  | neighbor-balance filter-pairs $chromsizes --regions $captureregions \
  | gzip > ${name}.nodups.select.shifted.pairs.gz"


####################################################################################################
# Load into coolers and balance.
sbatch --mem=32G --cpus-per-task 16 -J $name --dependency=SINGLETON -p $partition \
  --wrap="zcat ${name}.nodups.select.shifted.pairs.gz \
  | cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 $chromsizes:200 - ${name}.cool & \
  zcat ${name}.nodups.select.shifted.pairs.gz \
  | neighbor-balance filter-pairs $chromsizes --direction inward \
  | cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 $chromsizes:200 - ${name}_inward.cool & \
  wait;"

sbatch --mem=32G --cpus-per-task 16 -J $name --dependency=SINGLETON -p $partition \
  --wrap="neighbor-balance remove-inward ${name}.cool ${name}_inward.cool ${name}_minus_inward.cool"

sbatch --mem=32G --cpus-per-task 16 -J $name --dependency=SINGLETON -p $partition \
  --wrap="cooler zoomify --balance --balance-args '--ignore-diags 0' --resolutions 200,400,800,1600,3200,6400,12800,25600,51200,102400 \
          --out ${name}_minus_inward.mcool ${name}_minus_inward.cool"


#####################################################################################################
# Add neighbor balance weights to cooler
# This is useful for genome-wide datasets, for region-capture datasets process regions individually (see README.md).

sbatch --mem=32G -J $name --dependency=SINGLETON -p $partition \
  --wrap="neighbor-balance neighbor-balance-cooler ${name}_minus_inward.mcool"


####################################################################################################
# Merge replicates.
# Use the below example command to produce a merged cool file and then use the above commands to produce balanced files.
name=WT
sbatch -J name --wrap="cooler merge ${name}_minus_inward.cool \
  WT_BR1/select_corrected_minus_inward.mcool::/resolutions/200 \
  WT_BR2/WT_BR2_minus_inward.mcool::/resolutions/200"
