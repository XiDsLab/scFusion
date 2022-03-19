from __future__ import print_function
import sys
import codecs


# ***** readme *****
# This code remove the low mappability records of Chimeric sam by star
# if the mappability < 1, delete the record locating at this position

ChimericSamFile = codecs.open(sys.argv[1], 'r',  encoding='utf-8', errors='ignore')
OutputFile = open(sys.argv[2], 'w')
chrflag = []
lastname = ''
linestore = []
flag = 0
aaaaa = 0
linecount = 0
lines = ChimericSamFile.readlines()
for line in lines:
    info = line.split('\t')
    if len(info) < 5:
        continue
    if len(line) > 10:
        if line[0] == '@':
            continue
        if info[2].find('M') > -1:
            continue
    if info[0] != lastname:
        bad = 0
        poss = []
        for item in linestore:
            info = item.split('\t')
            '''
            if info[5].find('N') > -1:
                bad = 1
                break
            '''
            poss.append(info[3])
        if len(poss) <= 1 or (len(poss) == 2 and abs(int(poss[0])-int(poss[1])) <= 10):
            continue
        for subline in linestore:
            try:
                sep = '\t'
                subinfo = subline.split('\t')
                if not subinfo[2].startswith('chr'):
                    subinfo[2] = 'chr' + subinfo[2]
                OutputFile.write(sep.join(subinfo))
            except:
                print(subline)
        info = line.split('\t')
        lastname = info[0]
        linestore = [line]
    else:
        linestore.append(line)
    if line == '':
        break
ChimericSamFile.close()
OutputFile.close()
