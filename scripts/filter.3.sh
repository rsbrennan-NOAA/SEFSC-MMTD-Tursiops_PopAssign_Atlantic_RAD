#!/bin/bash
#SBATCH --job-name=filter2
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-NC-PopulationAssignment-RAD/logout
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH --partition=standard
#SBATCH --time=40:00

source ~/.bashrc

mamba activate vcflib-1.0.9
module load bio/vcftools/0.1.16
module load bio/bcftools/1.11

INDIR=~/Tursiops-NC-PopulationAssignment-RAD/analysis/variants

vcftools --gzvcf ${INDIR}/filtered.6.vcf.gz --max-meanDP 139 --exclude-positions ~/Tursiops-NC-PopulationAssignment-RAD/scripts/HD_exclude.txt --recode --recode-INFO-all --stdout |  bgzip > ${INDIR}/filtered.final.vcf.gz

bcftools annotate --set-id '%CHROM\_%POS'  ${INDIR}/filtered.final.vcf.gz  -Oz -o ${INDIR}/filtered.final_ids.vcf.gz

bcftools index -t ${INDIR}/filtered.final.vcf.gz
bcftools index -t ${INDIR}/filtered.final_ids.vcf.gz




