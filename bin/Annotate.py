from __future__ import print_function
import sys
import pyensembl


myfile = open(sys.argv[1])
gtffile = sys.argv[2]
outfile = open(sys.argv[1][:-4] + '_geneanno.sam', 'w')
data = pyensembl.Genome(reference_name='GRCH37', annotation_name='my_genome_features', gtf_path_or_url=gtffile)
data.index()
with myfile:
    for line in myfile:
        if line[0] == '@':
            continue
        info = line.split('\t')
        chr = info[2]
        pos = int(info[3])
        nameresult = data.gene_names_at_locus(contig=chr, position=pos)
        k = 20
        if nameresult == []:
            nameresult = data.gene_names_at_locus(contig=chr, position=pos+k)
        if nameresult == []:
            nameresult = data.gene_names_at_locus(contig=chr, position=pos - k)
        genename = ''
        for item in nameresult:
            genename += item + ';'
        genename = genename[:-1]
        outfile.write(genename + '\t' + line)
