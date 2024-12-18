#!/bin/bash
#SBATCH --job-name=filter1
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-NC-PopulationAssignment-RAD/logout
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH --partition=standard
#SBATCH --time=4:00:00

source ~/.bashrc

mamba activate vcflib-1.0.9
module load bio/vcftools/0.1.16

INDIR=~/Tursiops-NC-PopulationAssignment-RAD/analysis/variants

echo "starting first missingness filter"

vcftools --gzvcf ${INDIR}/variants_raw_merged.vcf.gz \
	--max-missing 0.5 --recode-INFO-all --stdout | \
	bgzip > ${INDIR}/filtered.1.vcf.gz

NUM_VARIANTS1=$(zcat ${INDIR}/filtered.1.vcf.gz | grep -v '^#' | wc -l)
echo "Number of variants after missingness filter: ${NUM_VARIANTS1}"

vcfallelicprimitives ${INDIR}/variants_raw_merged.vcf.gz --keep-info --keep-geno | \
	vcfstreamsort |  bgzip > ${INDIR}/filtered.1.vcf.gz

NUM_VARIANTS1=$(zcat ${INDIR}/filtered.1.vcf.gz | grep -v '^#' | wc -l)
echo "Number of variants after 1st filter: ${NUM_VARIANTS1}"

echo "running vcfallelicprimitives"
vcfallelicprimitives ${INDIR}/filtered.1.vcf.gz --keep-info --keep-geno | \
	        vcfstreamsort |  bgzip > ${INDIR}/filtered.2.vcf.gz

NUM_VARIANTS2=$(zcat ${INDIR}/filtered.2.vcf.gz | grep -v '^#' | wc -l)
echo "Number of variants after primitives step: ${NUM_VARIANTS2}"



vcftools --gzvcf ${INDIR}/filtered.2.vcf.gz \
	        --mac 3 --remove-indels --max-alleles 2 --min-alleles 2 --minQ 30  \
		--recode-INFO-all --stdout | \
		vcftools --vcf - --minDP 3 --recode-INFO-all --stdout | \
		bgzip > ${INDIR}/filtered.3.vcf.gz

NUM_VARIANTS3=$(zcat ${INDIR}/filtered.3.vcf.gz | grep -v '^#' | wc -l)
echo "Number of variants after quality filters: ${NUM_VARIANTS3}"

# stricter missingness
vcftools --gzvcf ${INDIR}/filtered.3.vcf.gz --max-missing 0.6 --recode-INFO-all --stdout | \
	    bgzip > ${INDIR}/filtered.4.vcf.gz

NUM_VARIANTS4=$(zcat ${INDIR}/filtered.4.vcf.gz | grep -v '^#' | wc -l)
echo "Number of variants after stricter missingness filter: ${NUM_VARIANTS4}"

# check missingness:
cd ${INDIR}
vcftools --gzvcf ${INDIR}/filtered.4.vcf.gz --missing-indv


