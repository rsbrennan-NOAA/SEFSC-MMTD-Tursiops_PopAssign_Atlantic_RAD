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
#SBATCH --time=7:00:00

source ~/.bashrc

mamba activate vcflib-1.0.9
module load bio/vcftools/0.1.16
module load bio/bcftools/1.11
module load bio/htslib/1.19

INDIR=~/Tursiops-NC-PopulationAssignment-RAD/analysis/variants

echo "starting first missingness filter"

#bcftools filter -S . -e 'FMT/DP<3' ${INDIR}/variants_raw_merged.vcf.gz | \
#	bcftools view -i 'F_MISSING < 0.4' -Oz -o ${INDIR}/filtered.1.vcf.gz
bcftools view -i 'INFO/DP/N_SAMPLES >= 3 && F_MISSING < 0.4' -Oz -o ${INDIR}/filtered.1.vcf.gz ${INDIR}/variants_raw_merged.vcf.gz

# Index the output
#bcftools index -t ${INDIR}/filtered.1.vcf.gz

#vcftools --gzvcf ${INDIR}/variants_raw_merged.vcf.gz \
#	--minDP 3 --recode --recode-INFO-all --stdout | \
#	vcftools --vcf - --max-missing 0.6 --recode --recode-INFO-all --stdout | \
#	bgzip > ${INDIR}/filtered.1.vcf.gz

tabix -p vcf -f ${INDIR}/filtered.1.vcf.gz

echo "done with first missingness filter"

NUM_VARIANTS1=$(zcat ${INDIR}/filtered.1.vcf.gz | grep -v '^#' | wc -l)
echo "Number of variants after first missingness filter: ${NUM_VARIANTS1}"

echo "running vcfallelicprimitives"
vcfallelicprimitives ${INDIR}/filtered.1.vcf.gz --keep-info --keep-geno | \
	        vcfstreamsort |  bgzip > ${INDIR}/filtered.2.vcf.gz

echo "done with allelic primitive"
NUM_VARIANTS2=$(zcat ${INDIR}/filtered.2.vcf.gz | grep -v '^#' | wc -l)
echo "Number of variants after primitives step: ${NUM_VARIANTS2}"

echo "start 2nd filtering"
vcftools --gzvcf ${INDIR}/filtered.2.vcf.gz \
	        --mac 3 --remove-indels --max-alleles 2 --min-alleles 2 --minQ 30  \
		--recode --recode-INFO-all --stdout | \
		bgzip > ${INDIR}/filtered.3.vcf.gz

echo "finish second filtering"

NUM_VARIANTS3=$(zcat ${INDIR}/filtered.3.vcf.gz | grep -v '^#' | wc -l)
echo "Number of variants after second filters: ${NUM_VARIANTS3}"


# check missingness:
cd ${INDIR}
vcftools --gzvcf ${INDIR}/filtered.3.vcf.gz --missing-indv


