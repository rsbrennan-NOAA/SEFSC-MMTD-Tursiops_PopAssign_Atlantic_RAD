#!/bin/bash
#SBATCH --job-name=pca
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-NC-PopulationAssignment-RAD/logout
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH --partition=standard
#SBATCH --time=1:00:00

source ~/.bashrc

module load bio/plink/1.90b6.23
module load bio/bcftools/1.11
module load bio/vcftools/0.1.16
module load bio/htslib/1.19

INDIR=~/Tursiops-NC-PopulationAssignment-RAD/analysis/

# add snp ids
bcftools annotate --set-id '%CHROM\_%POS' ${INDIR}/variants/filtered.final.vcf.gz -Oz -o ${INDIR}/variants/filtered.final_ids.vcf.gz

#make and move to structure directory
mkdir -p ${INDIR}/structure
cd  ${INDIR}/structure

plink --vcf ${INDIR}/variants/filtered.final_ids.vcf.gz \
	--indep-pairwise 50 5 0.2 --allow-extra-chr --double-id \
	--out variants_pruned

# make thinned vcf
vcftools --gzvcf ${INDIR}/variants/filtered.final_ids.vcf.gz --snps variants_pruned.prune.in  --recode --recode-INFO-all --stdout | \
       	bgzip > ${INDIR}/variants/filtered.final_LDthin.vcf.gz

# make plink files
plink --vcf ${INDIR}/variants/filtered.final_LDthin.vcf.gz \
	--extract variants_pruned.prune.in \
	--make-bed --out variants_NoLD \
	--allow-extra-chr --double-id

plink --bfile variants_NoLD \
	--pca --out all_indivs_PCA --allow-extra-chr --double-id
i
