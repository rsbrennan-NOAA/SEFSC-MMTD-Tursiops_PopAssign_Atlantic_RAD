#!/bin/bash
#SBATCH -D /home/rbrennan/spermWhaleRad/logout
#SBATCH --mail-type=END
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --partition=standard
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=80G
#SBATCH --time=05:00:00
#SBATCH --job-name=bwa_index
#SBATCH --output=%x.%A.%a.out
#SBATCH --error=%x.%A.%a.err

module load aligners/bwa-mem2/2.2.1

cd ~/reference_genomes/bottlenose_dolphin/

bwa-mem2 index ~/reference_genomes/bottlenose_dolphin/GCF_011762595.1_mTurTru1.mat.Y_genomic.fna.gz
