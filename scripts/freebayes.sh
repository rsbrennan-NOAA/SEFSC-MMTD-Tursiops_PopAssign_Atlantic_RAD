#!/bin/bash 
#SBATCH --job-name=freebayes_parallel
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-NC-PopulationAssignment-RAD/logout
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --partition=standard
#SBATCH --time=30:00:00
#SBATCH --array=1-362%8 

##################################
# run freebayes to call variants
##################################

#load software---------------------------------------------------------------
source ~/.bashrc

mamba activate freebayes-1.3.6 

#input, output files, directories--------------------------------------------
INDIR=~/Tursiops-NC-PopulationAssignment-RAD/analysis/merged_bams
OUTDIR=~/Tursiops-NC-PopulationAssignment-RAD/analysis/variants/chromosomes

#reference genome
GEN=~/reference_genomes/bottlenose_dolphin/GCF_011762595.1_mTurTru1.mat.Y_genomic.fna

# input files
BAMLIST=~/Tursiops-NC-PopulationAssignment-RAD/scripts/bam.list

# Get chromosome name from the array task ID
CHROM=$(awk "NR==${SLURM_ARRAY_TASK_ID}" ${GEN}.fai | cut -f1)

# these are freebayes scripts found in the same location as the executable
# Make regions file for this chromosome only
MAKEREGIONS=~/bin/freebayes-1.3.8/scripts/fasta_generate_regions.py
REG=$OUTDIR/regions_${CHROM}.txt
python3 $MAKEREGIONS ${GEN}.fai 5000000 | grep "^${CHROM}" > $REG

# Run freebayes-parallel for only each chromosome

fb_parallel=~/bin/freebayes-1.3.8/scripts/freebayes-parallel

bash $fb_parallel \
	$REG 20 \
	-f ${GEN} \
	--bam-list $BAMLIST \
	-m 30 \
	-q 20 \
    --use-best-n-alleles 4 \
    --min-alternate-count 2 | \
bgzip -c > $OUTDIR/variants_raw_${CHROM}.vcf.gz

# index the vcf files
tabix -p vcf $OUTDIR/variants_raw_${CHROM}.vcf.gz
