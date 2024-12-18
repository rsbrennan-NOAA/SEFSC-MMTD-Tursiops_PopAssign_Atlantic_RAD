#!/bin/bash
#SBATCH --job-name=vcf_concat
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-NC-PopulationAssignment-RAD/logout
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --partition=standard
#SBATCH --time=4:00:00

# load modules
module load bio/bcftools/1.11
module load bio/htslib/1.19
OUTDIR=~/Tursiops-NC-PopulationAssignment-RAD/analysis/variants
INDIR=~/Tursiops-NC-PopulationAssignment-RAD/analysis/variants/chromosomes
# make list of VCF files
ls $INDIR/variants_raw_*.vcf.gz > $OUTDIR/vcf_list.txt

# merge vcfs
bcftools concat -f $OUTDIR/vcf_list.txt -O z -o $OUTDIR/variants_raw_merged.vcf.gz

# index
tabix -p vcf $OUTDIR/variants_raw_merged.vcf.gz
