#!/usr/bin/env python

import gzip
import bioinfo
import matplotlib.pyplot as plt
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="Reads FASTA file and concatenate all sequence lines into one single line")
    parser.add_argument("-f", "--filename", help="Path to file", required=True)
    parser.add_argument("-o", "--outputfile", help="Path to output file", type=str)
    parser.add_argument("-u", "--outputimage", help="Path to output histogram image", type=str)
    return parser.parse_args()
args=get_args()	

f = args.filename  # input file name
o = args.outputfile # output file name
u = args.outputimage 



# Create an empty list:
qscore_sum = []
i = 0
while i <= 100:
    qscore_sum.append(0.0)
    i+=1

# Using empty list to calculate a running sum of qscores for each base position:
# with gzip.open (f, "r") as file:
with open (f,"r") as file:
    i = 0
    num_lines = 0
    for line in file:
        line = line.strip()
        i+=1
        num_lines+=1
        if i%4 == 0:
            for base in range (len(line)):
                qscore_sum [base] += bioinfo.convert_phred(line[base])

# Use list of sums to calculate all the means:
num_record = num_lines/4
qscore_mean = []
for qscore in qscore_sum:
    qscore_mean.append(qscore/num_record)

# Create histogram and tsv file:
x = []
for i in range(len(qscore_mean)):
    x.append(i)
y = qscore_mean

with open (o,"w") as o:
    o.write("# Base\tMean Qscore\n")
    for i in range(len(qscore_mean)):
        o.write(f'{i}\t{qscore_mean[i]}\n')

plt.bar (x,y)
plt.xlabel("# of Base")
plt.ylabel("Mean QScore")
plt.title("Distribution of Quality Scores")
plt.savefig(u)



