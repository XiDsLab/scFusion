from __future__ import print_function
from __future__ import division
import sys
import random
import os
import pysam

# ***** readme *****
# This code extracts chimeric read from sam file for training, with pos and direction
# The input is *.sam

def ReverseComplement(str):
    return str[::-1].replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g').upper()


random.seed(1122)
chimericfile = open(sys.argv[1])
mappingpath = sys.argv[2]
linenum = len(chimericfile.readlines())
chimericfile.close()
count = 0
cellindex = []
for dir in os.listdir(mappingpath):
    cellindex.append(dir)
while count < linenum:
    thisindex = random.sample(cellindex, 1)[0]
    try:
        found = False
        for file in os.listdir(mappingpath + thisindex):
            if file.find('Aligned.sortedByCoord.out.bam') > -1:
                samfile = pysam.AlignmentFile(mappingpath + thisindex + '/' + file, 'rb')
                found = True
                break
    except:
        continue
    if not found:
        continue
    sam = []
    for r in samfile:
        sam.append(r)
    thiscount = 0
    while thiscount < len(sam) / 5 and count < linenum:
        while True:
            a = random.randint(1, len(sam)) - 1
            b = random.randint(1, len(sam)) - 1
            chr1 = str(sam[a].reference_name)
            chr2 = str(sam[b].reference_name)
            allread1 = sam[a].seq
            allread2 = sam[b].seq
            readlength = len(allread1)
            if not chr1.startswith('chr'):
                chr1 = 'chr' + chr1
            if not chr2.startswith('chr'):
                chr2 = 'chr' + chr2
            if len(sam[a].cigar) > 1 or len(sam[b].cigar) > 1:
                continue
            if not (chr1 == '*' or chr2 == '*' or chr1 == 'chrM' or chr2 == 'chrM'):
                break
        read1length = 30
        read2length = 60 - read1length
        c = random.randint(0, 60 - read1length - 1)
        d = random.randint(0, 60 - read2length - 1)
        try:
            read1 = sam[a].seq[c:c + read1length]
        except:
            sys.stderr.write(str(sam[a]))
        line2 = sam[b]
        try:
            read2 = sam[b].seq[d:d + read2length]
        except:
            sys.stderr.write(str(sam[b]))
        e = random.randint(0, 1)
        f = random.randint(0, 1)
        if e == 0:
            e = -1
            read2 = ReverseComplement(read2)
            pos2 = sam[b].pos + d + read2length - 1
        else:
            pos2 = sam[b].pos + d
        if f == 0:
            f = -1
            read1 = ReverseComplement(read1)
            pos1 = sam[a].pos + c
        else:
            pos1 = sam[a].pos + c + read1length - 1


        if f == -1:
            direct1 = '+'
        else:
            direct1 = '-'
        if e == 1:
            direct2 = '+'
        else:
            direct2 = '-'
        if read1.find('N') == -1 and read2.find('N') == -1 and len(read1 + read2) == 60:
            print(read1.upper() + read2.upper() + '\t' + str(read1length) + '\t', end='')
            print(chr1 + ':' + str(pos1) + ':' + direct1 + '\t' + chr2 + ':' + str(pos2) + ':' + direct2)
            count += 1
            thiscount += 1
    samfile.close()