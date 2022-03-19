from __future__ import print_function
from __future__ import division
import sys
import math
import random


# ***** readme *****
# This code extracts chimeric read from sam file for training, with pos and direction
# The input is *ChiDist_middle.txt
random.seed(1122)
infile = open(sys.argv[1])
lines = infile.readlines()
totallines = len(lines)
uselines = random.sample(lines[int(totallines/10):], math.floor(min(int(totallines * 0.4), max(15000, totallines/9))))
for line in uselines:
    info = line.rstrip().split('\t')
    gene1 = info[0]
    gene2 = info[1]
    if gene1.startswith('IG')  or gene2.startswith('IG')  or gene2.startswith('TRA') or gene1.startswith('TRA'):
        continue
    print(info[-4] + '\t' + info[-3] + '\t' + info[6] + ':' + info[8] + ':' + info[-2] + '\t' + info[7] + ':' + info[9] + ':' + info[-1])
infile.close()
