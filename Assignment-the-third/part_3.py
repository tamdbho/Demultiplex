#!/usr/bin/env python
import argparse
import bioinfo
import gzip
import re


def get_args():
    parser = argparse.ArgumentParser(description="Reads FASTA file and concatenate all sequence lines into one single line")
    parser.add_argument("-r1", "--read1file", help="Path to biological read file 1 - this would be R1", type=str)
    parser.add_argument("-r2", "--read2file", help="Path to biological read file 2 - this would be R4", type=str)
    parser.add_argument("-i1", "--index1file", help="Path to index file 1 - this would be R2", type=str)
    parser.add_argument("-i2", "--index2file", help="Path to index file 2 - this would be R3", type=str)
    parser.add_argument("-e", "--expectedindex", help="Path to file containing expected indexes", type=str)
    parser.add_argument("-o", "--outputdir", help="Path to where the output files will be located", type=str)
    parser.add_argument("-c", "--qualityscore", help="Quality score cutoff", type=int)
    return parser.parse_args()

args=get_args()	

read1_file = args.read1file
read2_file = args.read2file
index1_file = args.index1file
index2_file = args.index2file
expected_index = args.expectedindex
cutoff = args.qualityscore
output = args.outputdir

# Initiate some counter variable for a report file: count number of unknown reads, index-hopped reads, matched reads per index pair
unknown_count = 0 
index_hop_count = 0
dual_matched_count = 0
matched_count = {}  # Create dictionary to report count and percentage of each dual-matched pair
num_lines = 0
ihop_pair_count = {}    # Create dictionary to report count and percentage of each index-hopped pair


def record_perloop (line) -> tuple:
    '''To be used while iterating over FASTQ file to temporary store each component of a record in a tuple 
    - do not do the looping - only saving variables until the next loop'''
    header = line.readline().strip()
    sequence = line.readline().strip()
    comment = line.readline().strip()
    quality = line.readline().strip()
    record = (header,sequence,comment,quality)
    return record

# Create a set containing all known index using the indexes.txt file:
known_index = set()
with open(expected_index,"r") as file:
    line = file.readline().strip()
    for line in file:
        line = line.split()
        known_index.add(line[-1])

# Create a function that will check quality score of the indexes - to be used in one of the conditions:
# def mean_qscore (score):
#     mean_q = []
#     for q in score:
#         mean_q.append(bioinfo.convert_phred(q))
#     mean = sum(mean_q)/len(mean_q)
#     return mean

# Intitialize a dictionary to keep track of index pair (key) and their count (value) - will need this number for analysis later
write_dict = {}

# Initialize all files needed for writing from the start using open() so it does not autimatically closes. 
indexhop_R1 = open(f'{output}Index_hopped.R1.fq.out',"w")
indexhop_R2 = open(f'{output}Index_hopped.R2.fq.out',"w")
write_dict["index_hop"] = [indexhop_R1, indexhop_R2]

unknown_R1 = open(f'{output}Unknown.R1.fq.out',"w")
unknown_R2 = open(f'{output}Unknown.R2.fq.out',"w")
write_dict["unknown"] = [unknown_R1,unknown_R2]

i=0
for index in known_index:
    i+=1
    fh1 = open(f'{output}{index}-{index}.R1.fq.out',"w") 
    fh2 = open(f'{output}{index}-{index}.R2.fq.out',"w")
    write_dict[index] = [fh1,fh2]
    # loop through set of index, create a file to write to for each index
    # at the same time: append the index as key to dictionary and file handles as values so we can call on the file later in our loop to write

with gzip.open (read1_file,"rt") as r1, gzip.open(read2_file,"rt") as r2, gzip.open(index1_file,"rt") as i1, gzip.open(index2_file,"rt") as i2:
# with open (read1_file,"rt") as r1, open(read2_file,"rt") as r2, open(index1_file,"rt") as i1, open(index2_file,"rt") as i2:
    i = 0
    while True:
        record_r1 = record_perloop(r1)
        record_r2 = record_perloop(r2)
        record_i1 = record_perloop(i1)
        record_i2 = record_perloop(i2)
        # each record is a tuple (header,sequence,comment,quality score)
        if record_r1 == ("","","",""):
            break
        i+=1
        num_lines += 1
        rc_i2 = bioinfo.rev_comp_DNA(record_i2[1])
        if record_i1[1] not in known_index or rc_i2 not in known_index:
            write_dict["unknown"][0].write(f'{record_r1[0]} - {record_i1[1]}-{rc_i2}\n{record_r1[1]}\n{record_r1[2]}\n{record_r1[3]}\n')
            write_dict["unknown"][1].write(f'{record_r2[0]} - {record_i1[1]}-{rc_i2}\n{record_r2[1]}\n{record_r2[2]}\n{record_r2[3]}\n')
            unknown_count += 1
        else:
            if record_i1[1] == rc_i2:
                write_dict[record_i1[1]][0].write(f'{record_r1[0]} - {record_i1[1]}-{rc_i2}\n{record_r1[1]}\n{record_r1[2]}\n{record_r1[3]}\n')
                write_dict[record_i1[1]][1].write(f'{record_r2[0]} - {record_i1[1]}-{rc_i2}\n{record_r2[1]}\n{record_r2[2]}\n{record_r2[3]}\n')
                if not matched_count.get(record_i1[1]):
                    matched_count[record_i1[1]] = 1
                else:
                    matched_count[record_i1[1]] += 1
                dual_matched_count += 1
            else:
                write_dict["index_hop"][0].write(f'{record_r1[0]} - {record_i1[1]}-{rc_i2}\n{record_r1[1]}\n{record_r1[2]}\n{record_r1[3]}\n')
                write_dict["index_hop"][1].write(f'{record_r2[0]} - {record_i1[1]}-{rc_i2}\n{record_r2[1]}\n{record_r2[2]}\n{record_r2[3]}\n')
                ihop_pair = record_i1[1] + " - " + rc_i2
                if not ihop_pair_count.get(ihop_pair):
                    ihop_pair_count[ihop_pair] = 1
                else:
                    ihop_pair_count[ihop_pair] += 1
                index_hop_count += 1

# Create a report file:
percent_matched = dual_matched_count/num_lines*100
percent_hop = index_hop_count/num_lines*100
percent_unknown = unknown_count/num_lines*100

with open (f'{output}Final_Report.tsv',"w") as report:
    report.write("Input files: \n")
    report.write(f'{read1_file}\n{read2_file}\n{index1_file}\n{index2_file}\n\n')
    
    report.write("SUMMARY\n")
    report.write (f'Total reads = {num_lines}\n')
    report.write(f'Total dual-matched count = {dual_matched_count}\tPercent = {percent_matched}\n')
    report.write (f'Total unknown reads = {unknown_count}\tPercent = {percent_unknown}\n')
    report.write (f'Total index-hopped reads = {index_hop_count}\tPercent = {percent_hop}\n\n')
    
    report.write("Dual-matched pairs:\n")
    report.write("Index pair\tCount\tPercent\n")
    for index, count in matched_count.items():
        percent = count/num_lines*100
        report.write(f'{index}\t{count}\t{percent}\n')
    report.write("\n")
    
    report.write("Index-hopped pairs:\n")
    report.write("Index-hopped\tCount\tPercent\n")
    for ihop, counter in ihop_pair_count.items():
        ihop_pair_percent = counter/num_lines*100
        report.write(f'{ihop}\t{counter}\t{ihop_pair_percent}\n')

for index, fh in write_dict.items():
    fh[0].close()
    fh[1].close()