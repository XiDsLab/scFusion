from __future__ import print_function
from __future__ import division
import sys
import math
import os
import numpy

# ***** readme *****
# 这个代码将所有的fastq文件重命名为数字系列，并且保存映射表
# 也可以反过来操作


folder = sys.argv[1]

mapfilepath = ''

if len(sys.argv) > 2:
    mapfilepath = sys.argv[2]

mapdict = {}
if mapfilepath == '':
    if os.path.exists(folder + '/RenameList.txt'):
        print('The input files may have already renamed, since RenameList.txt is detected in ' + folder + '. To force rename files, please remove RenameList.txt')
        exit(1)
    count = 0
    for filename in os.listdir(folder):
        if filename.endswith('_1.fastq'):
            count += 1
            os.rename(folder + '/' + filename, folder + '/' + str(count) + '_1.fastq')
            mapdict[count] = filename[:-8]
            try:
                os.rename(folder + '/' + filename[:-7] + '2.fastq', folder + '/' + str(count) + '_2.fastq')
            except Exception as e:
                sys.stderr.write(str(e) + ' Rename the pair: ' + filename[:-7] + '_2.fastq failed.\n')
        if filename.endswith('_1.fastq.gz'):
            count += 1
            os.rename(folder + '/' + filename, folder + '/' + str(count) + '_1.fastq.gz')
            mapdict[count] = filename[:-11]
            try:
                os.rename(folder + '/' + filename[:-10] + '2.fastq.gz', folder + '/' + str(count) + '_2.fastq.gz')
            except Exception as e:
                sys.stderr.write(str(e) + ' Rename the pair: ' + filename[:-10] + '_2.fastq.gz failed.\n')
    outfile = open(folder + '/RenameList.txt', 'w')
    for key in mapdict:
        outfile.write(str(key) + '\t' + mapdict[key] + '\n')
    outfile.close()
else:
    if not os.path.exists(folder + '/RenameList.txt'):
        print('RenameList.txt cannot be found in ' + folder + ', please check it')
        exit(1)
    mapfile = open(mapfilepath)
    for line in mapfile.readlines():
        info = line.rstrip().split('\t')
        mapdict[int(info[0])] = info[1]
        for key in mapdict:
            if os.path.exists(folder + '/' + str(key) + '_1.fastq'):
                os.rename(folder + '/' + str(key) + '_1.fastq', folder + '/' + str(mapdict[key]) + '_1.fastq')
            if os.path.exists(folder + '/' + str(key) + '_2.fastq'):
                os.rename(folder + '/' + str(key) + '_2.fastq', folder + '/' + str(mapdict[key]) + '_2.fastq')
            if os.path.exists(folder + '/' + str(key) + '_1.fastq.gz'):
                os.rename(folder + '/' + str(key) + '_1.fastq.gz', folder + '/' + str(mapdict[key]) + '_1.fastq.gz')
            if os.path.exists(folder + '/' + str(key) + '_2.fastq.gz'):
                os.rename(folder + '/' + str(key) + '_2.fastq.gz', folder + '/' + str(mapdict[key]) + '_2.fastq.gz')