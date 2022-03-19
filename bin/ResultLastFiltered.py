from __future__ import print_function
from __future__ import division
import sys
import os.path

resultfile = open(sys.argv[1])
geneposfile = open(sys.argv[2])
outfilepath = sys.argv[3]
annot = open(sys.argv[4])
geneposfile2path = sys.argv[5]
lncfilter = sys.argv[6]
nasfilter = sys.argv[7]
outfile = open(outfilepath, 'w')
trimmedfile = open(outfilepath + '.fail', 'w')
totalfile = open(outfilepath + '.total', 'w')
genedic = {}
linecount = 0
exondic = {}
annotgene = {}
LINCoption = True
NoApprovoption = True
for line in annot.readlines():
    if line.startswith('#'):
        continue
    info = line.split('\t')
    start = int(info[3])
    end = int(info[4])
    geneinfo = info[8].split('; ')
    genename = ''
    for item in geneinfo:
        if item.startswith('gene_name'):
            genename = item[11:-1]
            break
    if genename == '':
        for item in geneinfo:
            if item.startswith('gene_id'):
                genename = item[9:-1]
    if info[2] == 'exon':
        if genename != '':
            if genename not in exondic:
                exondic[genename] = []
            exondic[genename].append([start, end])
    if info[2] == 'gene':
        annotgene[genename] = [[start, end]]
for key in annotgene:
    if key not in exondic:
        exondic[key] = annotgene[key]
for line in geneposfile.readlines():
    linecount += 1
    info = line.rstrip().split('\t')
    if info[0] not in genedic:
        genedic[info[0]] = [linecount, info[1], int(info[2]), int(info[3]), info[4]]
if geneposfile2path != '':
    geneposfile2 = open(geneposfile2path)
    for line in geneposfile2.readlines():
        info = line.rstrip().split('\t')
        if info[0] not in genedic:
            genedic[info[0]] = [linecount, info[1], int(info[2]), int(info[3]), info[4]]
        else:
            if genedic[info[0]][4] != 'pseudogene' and genedic[info[0]][4] != 'lincRNA':
                genedic[info[0]][4] = info[4]
posdict = {}
for line in resultfile.readlines():
    if line.startswith('#'):
        continue
    info = line.split('\t')
    genes = info[0].split('--')
    suppcells = int(info[1])
    splitread = int(info[2])
    discordant = int(info[3])
    pos1 = int(info[4].split(':')[1])
    pos2 = int(info[5].split(':')[1])
    posdict[pos1 + pos2] = 1
    flag = ['\t', '\t', '\t', '\t', '\t', '\t']
    if NoApprovoption and ((genes[0].find('RP') > -1 and genes[0].find('-') > -1 and genes[0].find('.') > -1) or \
            (genes[1].find('RP') > -1 and genes[1].find('-') > -1 and genes[1].find('.') > -1) or \
            (genes[0].find('A') > -1 and genes[0].find('.') > -1) or \
            (genes[1].find('A') > -1 and genes[1].find('.') > -1) or \
            genes[0].find('MT-') > -1 or genes[1].find('MT-') > -1 or \
            (genes[0].find('CT') > -1 and genes[0].find('-') > -1 and genes[0].find('.') > -1) or \
            (genes[1].find('CT') > -1 and genes[1].find('-') > -1 and genes[1].find('.') > -1)):
        flag[0] = '\tRP'
    if discordant / splitread > 10 or suppcells >= 100 and splitread / suppcells < 1.6:
        flag[1] = '\tTooManyDiscordant'
    if genedic[genes[0]][1] == genedic[genes[1]][1] and splitread / suppcells < 5 and (abs(genedic[genes[0]][0] - genedic[genes[1]][0]) < 2 or not (genedic[genes[0]][2] > genedic[genes[1]][3] or genedic[genes[1]][2] > genedic[genes[0]][3])) or genes[0] == genes[1]:
        flag[2] = '\tOverlap'
    if genedic[genes[0]][1] == genedic[genes[1]][1] and abs(pos1 - pos2) < 200000 and splitread < 150 and info[0].find('IGH') == -1:
        flag[2] = '\tOverlap'
    if LINCoption and (genedic[genes[0]][4] == 'lincRNA' or genedic[genes[1]][4] == 'lincRNA' or genedic[genes[0]][4] == 'lncRNA' or \
            genedic[genes[1]][4] == 'lncRNA' or genes[0].startswith('LINC') or genes[1].startswith('LINC')):
        flag[3] = '\tlncRNA'
    if genedic[genes[0]][4] == 'pseudogene' or genedic[genes[1]][4] == 'pseudogene':
        flag[4] = '\tpseudogene'
    gene1exon = False
    gene2exon = False
    for item in exondic[genes[0]]:
        if pos1 - 100 <= item[1] and item[0] <= pos1 + 100:
            gene1exon = True
            break
    for item in exondic[genes[1]]:
        if pos2 - 100 <= item[1] and item[0] <= pos2 + 100:
            gene2exon = True
            break
    if not gene1exon or not gene2exon:
        flag[5] = '\tintron'
    if (flag[0] == '\t' or nasfilter == '1') and flag[1] == '\t' and flag[2] == '\t' and (flag[3] == '\t' or lncfilter == '1') and flag[4] == '\t' and flag[5] == '\t':
        outfile.write(line)
    else:
        trimmedfile.write(line.rstrip() + flag[0] + flag[1] + flag[2] + flag[3] + flag[4] + flag[5] + '\n')
    totalfile.write(line.rstrip() + flag[0] + flag[1] + flag[2] + flag[3] + flag[4]  + flag[5] + '\n')
