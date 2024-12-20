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
#SBATCH --time=3:00:00

source ~/.bashrc

mamba activate vcflib-1.0.9
module load bio/vcftools/0.1.16

INDIR=~/Tursiops-NC-PopulationAssignment-RAD/analysis/variants

vcftools --gzvcf ${INDIR}/filtered.3.vcf.gz --remove ~/Tursiops-NC-PopulationAssignment-RAD/scripts/rm_missing.txt  --recode --recode-INFO-all --stdout |  bgzip > ${INDIR}/filtered.4.vcf.gz 

vcftools --gzvcf ${INDIR}/filtered.4.vcf.gz  --max-missing 0.7 --maf 0.05 --min-meanDP 10 --recode --recode-INFO-all  --stdout >  ${INDIR}/filtered.5.vcf

vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" ${INDIR}/filtered.5.vcf |bgzip > ${INDIR}/filtered.6.vcf.gz

#depths
vcftools --gzvcf ${INDIR}/filtered.6.vcf.gz --site-mean-depth --out ${INDIR}/filtered.6.depth
