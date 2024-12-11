#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-NC-PopulationAssignment-RAD/logout
#SBATCH -c 14
#SBATCH --partition=standard
#SBATCH --time=6:00:00

source ~/.bashrc

# load fastqc and multiqc
mamba activate multiqc-1.17

module load bio/fastqc/0.11.9


cd ~/Tursiops-NC-PopulationAssignment-RAD/

for input_dir in 1842 3420 3421 3422 3423 3424 3425
do
    mkdir -p analysis/fastqc_rawData/${input_dir}/
    fastqc -t 14 --noextract -o analysis/fastqc_rawData/${input_dir}/ data/${input_dir}/*.fastq.gz

    multiqc analysis/fastqc_rawData/${input_dir}/ -o  analysis/fastqc_rawData/ --filename multiqc_report_${input_dir}
done

