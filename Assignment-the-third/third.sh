#!/bin/bash
#SBATCH --account=bgmp                    
#SBATCH --partition=compute               
#SBATCH --mail-type=ALL                  
#SBATCH --cpus-per-task=8                 
#SBATCH --mem=32GB  

conda activate bgmp_py311
read_file1=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz
read_file2=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz

index_file1=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
index_file2=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz

known_index=/projects/bgmp/tamho/bioinfo/Bi622/Demultiplex/indexes.txt

output_dir=/projects/bgmp/tamho/bioinfo/Bi622/Demultiplex/Assignment-the-third/test_output_noqc_new/



/usr/bin/time -v ./part_3.py -r1 $read_file1 -r2 $read_file2 -i1 $index_file1 -i2 $index_file2 -e $known_index -o $output_dir 
