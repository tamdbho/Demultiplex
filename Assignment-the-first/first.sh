#!/bin/bash
#SBATCH --account=bgmp                    
#SBATCH --partition=compute               
#SBATCH --mail-type=ALL                  
#SBATCH --cpus-per-task=8                 
#SBATCH --mem=32GB  

conda activate bgmp_py311
file="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
tsv="Read1.tsv"
png="Read1.png"

/usr/bin/time -v ./first.py -f $file -o $tsv -u $png -l 101
