import os
import sys
import argparse
import subprocess

lastd = sys.argv[0].rfind('/')
if lastd == -1:
    programdir = './'
else:
    programdir = sys.argv[0][:lastd] + '/'
parser = argparse.ArgumentParser(description='Single-cell Gene Fusion Detection (c) Zjie Jin', add_help=False)
subparsers = parser.add_subparsers(help='Select the function', required=True)
parser1 = subparsers.add_parser('BuildSTARIndex', help='Build the STAR Index', add_help=False)
group11 = parser1.add_argument_group('Required parameters')
group12 = parser1.add_argument_group('Optional parameters')
group11.add_argument("-s", "--STARIndex", help='The STAR index output folder.')
group11.add_argument("-g", "--Genome", help='The reference file (*.fasta or *.fa).')
group11.add_argument("-a", "--Annotation", help='The gtf annotation file (*.gtf).')
group12.add_argument("-t", "--Thread", default='8', help='Number of threads can be used, default is 8.')
parser2 = subparsers.add_parser('Rename', help='Rename files using consecutive numbers', add_help=False)
group2 = parser2.add_argument_group('Required parameters')
group2.add_argument("-f", "--FileDir", help='The folder of input data.')
parser3 = subparsers.add_parser('RestoreName', help='Restore file names.', add_help=False)
group3 = parser3.add_argument_group('Required parameters')
group3.add_argument("-f", "--FileDir", help='The folder of input data.')
parser4 = subparsers.add_parser('Prepare', help='Generate some necessary files.', add_help=False)
group4 = parser4.add_argument_group('Required parameters')
group4.add_argument("-d", "--GenomeDir", help='The path of folder saving the generated files [GENOMEDIR].')
group4.add_argument("-g", "--Genome", help='The reference file (*.fasta or *.fa).')
group4.add_argument("-a", "--Annotation", help='The gtf annotation file (*.gtf).')
parser5 = subparsers.add_parser('ReadMapping', help='Mapping pair-end reads using STAR.', add_help=False)
group51 = parser5.add_argument_group('Required parameters')
group52 = parser5.add_argument_group('Optional parameters')
group51.add_argument("-f", "--FileDir", help='The folder of input data.')
group51.add_argument("-b", "--Begin", help='The first index of input files.')
group51.add_argument("-e", "--End", help='The last index of input files.')
group51.add_argument("-s", "--STARIndex", help='The STAR index folder.')
group51.add_argument("-o", "--OutDir", help='The output folder of the results and temporal files.')
group52.add_argument("-t", "--Thread", default='8', help='Number of threads can be used, default is 8')
parser6 = subparsers.add_parser('ReadProcessing', help='Process the chimeric reads.', add_help=False)
group61 = parser6.add_argument_group('Required parameters')
group62 = parser6.add_argument_group('Optional parameters')
group61.add_argument("-d", "--GenomeDir", help='The path of GENOMEDIR.')
group61.add_argument("-o", "--OutDir", help='The output folder of the results and temporal files.')
mappabilitygroup = group62.add_mutually_exclusive_group()
mappabilitygroup.add_argument("-m", "--Mappability", default=programdir + '/data/hg19mappability75.txt',
                              help='The mappability file, default is hg19mappability75.txt in the data folder.')
mappabilitygroup.add_argument("-M", "--NoMappabilityFilter", help='Keep all reads and do not apply mappability filter.',
                              action='store_true')
group61.add_argument("-b", "--Begin", help='The first index of input files.')
group61.add_argument("-e", "--End", help='The last index of input files.')
parser7 = subparsers.add_parser('FusionCandidate', help='Identify the fusion candidates', add_help=False)
group71 = parser7.add_argument_group('Required parameters')
group72 = parser7.add_argument_group('Optional parameters')
group71.add_argument("-o", "--OutDir", help='The output folder of the results and temporal files.')
group71.add_argument("-b", "--Begin", help='The first index of input files.')
group71.add_argument("-e", "--End", help='The last index of input files.')
group71.add_argument("-d", "--GenomeDir", help='The path of GENOMEDIR.')
group72.add_argument("-p", "--Prefix", default='.',
                     help='The prefix of result file, default is blank. This should be specified if users want to compare the results of different settings.')
parser8 = subparsers.add_parser('Retrain', help='Retrain the network using current data.', add_help=False)
group81 = parser8.add_argument_group('Required parameters')
group82 = parser8.add_argument_group('Optional parameters')
group81.add_argument("-o", "--OutDir", help='The output folder of the results and temporal files.')
group82.add_argument("-p", "--Prefix", default='.',
                     help='The prefix of result file, default is blank. This should be specified if users want to compare the results of different settings.')
group82.add_argument("-w", "--Weight", default=programdir + '/data/weight-V9-2.hdf5',
                     help='The initial weight file of the deep-learning network, default is "weight-V9-2.hdf5" in the data folder.')
group82.add_argument("-c", "--Epoch", default='10',
                     help='The number of epochs in the retraining step, default is 10')
parser9 = subparsers.add_parser('ArtifactScoring', help='Find the artifacts', add_help=False)
group91 = parser9.add_argument_group('Required parameters')
group92 = parser9.add_argument_group('Optional parameters')
group91.add_argument("-d", "--GenomeDir", help='The path of GENOMEDIR.')
group91.add_argument("-o", "--OutDir", help='The output folder of the results and temporal files.')
group91.add_argument("-b", "--Begin", help='The first index of input files.')
group91.add_argument("-e", "--End", help='The last index of input files.')
group92.add_argument("-p", "--Prefix", default='.',
                     help='The prefix of result file, default is blank. This should be specified if users want to compare the results of different settings.')
group92.add_argument("-w", "--Weight", default=programdir + '/data/weight-V9-2.hdf5',
                     help='The weight file used in the deep-learning network, default is "weight-V9-2.hdf5" in the data folder.')
parser10 = subparsers.add_parser('FusionReport', help='Report gene fusions.', add_help=False)
group101 = parser10.add_argument_group('Required parameters')
group102 = parser10.add_argument_group('Optional parameters')
group101.add_argument("-f", "--FileDir", help='The folder of input data.')
group101.add_argument("-d", "--GenomeDir", help='The path of GENOMEDIR.')
group101.add_argument("-o", "--OutDir", help='The output folder of the results and temporal files.')
group101.add_argument("-b", "--Begin", help='The first index of input files.')
group101.add_argument("-e", "--End", help='The last index of input files.')
group102.add_argument("-p", "--Prefix", default='.',
                      help='The prefix of result file, default is blank. This should be specified if users want to compare the results of different settings.')
group102.add_argument("-v", "--PvalueCutoff", default='0.05', help='Pvalue(FDR) cutoff of the statistical model, default is 0.05.')
group102.add_argument("-n", "--ArtifactScoreCutoff", default='0.75', help='Artifact score cutoff, default is 0.75')
group102.add_argument("--LncRNAFilterOff", help='turn off the LncRNA filter', action='store_true')
group102.add_argument("--NoApprovedSymbolFilterOff", help='turn off the no-approved-symbol filter', action='store_true')
parser11 = subparsers.add_parser('DeleteTempFiles', help='Delete temporary files generated by scFusion. ', add_help=False)
group111 = parser11.add_argument_group('Required parameters')
group112 = parser11.add_argument_group('Optional parameters')
group112.add_argument("--AllTempFiles", help='Delete all temporary files generated by scFusion and STAR.', action='store_true')
group112.add_argument("--AllUnimportantTempFiles", help='Delete all unimportant temporary files generated by scFusion.', action='store_true')
group112.add_argument("--STARMappingFiles", help='Delete STAR mapping results.', action='store_true')
group111.add_argument("-o", "--OutDir", help='The output folder of the results and temporal files.')


if len(sys.argv) == 1:
    parser.print_help()
    exit()

args = parser.parse_args()

if sys.argv[1] == 'Rename':
    if not args.FileDir:
        parser2.print_help()
        print('Please specify all required parameters!')
        exit()
    aa = os.system('python ' + programdir + '/bin/RenameFastqFiles.py ' + args.FileDir)
    if aa == 0:
        print('Successfully rename files to consecutive numbers!')
        print('Please keep RenameList.txt carefully, or the file names can not be restored.')
    else:
        print('ERROR!!!!!')

if sys.argv[1] == 'RestoreName':
    if not args.FileDir:
        parser3.print_help()
        print('Please specify all required parameters!')
        exit()
    aa = os.system(
        'python ' + programdir + '/bin/RenameFastqFiles.py ' + args.FileDir + ' ' + args.FileDir + '/RenameList.txt')
    if aa == 0:
        print('Successfully restore file names!')
    else:
        print('ERROR!!!!!')

if sys.argv[1] == 'BuildSTARIndex':
    if not args.STARIndex or not args.Genome or not args.Annotation:
        parser1.print_help()
        print('Please specify all required parameters!')
        exit()
    os.system('STAR --runMode genomeGenerate --genomeDir ' + args.STARIndex + ' --runThreadN ' +
              args.Thread + ' --genomeFastaFiles ' + args.Genome + ' --sjdbGTFfile ' +
              args.Annotation)

if sys.argv[1] == 'Index':
    if not args.GenomeDir or not args.Genome or not args.Annotation:
        parser4.print_help()
        print('Please specify all required parameters!')
        exit()
    aa = os.system('mkdir -p ' + args.GenomeDir + ' && cp ' + args.Genome + ' ' + args.GenomeDir + '/ref.fa && cp ' +
                   args.Annotation + ' ' + args.GenomeDir + '/ref_annot.gtf && python ' + programdir +
                   '/bin/Addchr2gtf.py ' + args.GenomeDir + '/ref_annot.gtf > ' + args.GenomeDir +
                   '/ref_annot.gtf.added && python ' + programdir + '/bin/GetExonPos.py ' + args.GenomeDir +
                   '/ref_annot.gtf.added > ' + args.GenomeDir + '/exon_probe.hg19.gene.new.bed && python ' + programdir +
                   '/bin/GetGenePos.py ' + args.GenomeDir + '/ref_annot.gtf.added > ' + args.GenomeDir +
                   '/GenePos.txt && pyensembl install --reference-name GRCH37 --annotation-name my_genome_features --gtf '
                   + args.GenomeDir + '/ref_annot.gtf.added')
    if aa == 0:
        print('Finish Indexing!')
    else:
        print('ERROR!!!!!')

if sys.argv[1] == 'ReadMapping':
    if not args.FileDir or not args.OutDir or not args.Begin or not args.End or not args.STARIndex:
        parser5.print_help()
        print('Please specify all required parameters!')
        exit()
    aa = os.system(
        'sh ' + programdir + '/bin/StarMapping_Chimeric.sh ' + args.FileDir + ' ' + args.Begin + ' ' + args.End + ' ' + args.OutDir +
        'STARMapping/ ' + args.STARIndex + ' ' + args.Thread)
    if aa == 0:
        print('Finish mapping! Index: ' + args.Begin + ' ~ ' + args.End)
    else:
        print('ERROR!!!!!')

if sys.argv[1] == 'ReadProcessing':
    if not args.OutDir or not args.Begin or not args.End or not args.GenomeDir or not args.NoMappabilityFilter and not args.Mappability:
        parser6.print_help()
        print('Please specify all required parameters!')
        exit()
    if args.NoMappabilityFilter:
        aa = os.system(
            'sh ' + programdir + '/bin/CombinePipeline_before_FS_NoMappa.sh ' + args.OutDir + ' ' + args.Begin + ' ' + args.End + ' ' +
            args.GenomeDir + '/ref_annot.gtf.added ' + args.GenomeDir + '/exon_probe.hg19.gene.new.bed ' + programdir + '/bin/')
    else:
        aa = os.system(
            'sh ' + programdir + '/bin/CombinePipeline_before_FS.sh ' + args.OutDir + ' ' + args.Begin + ' ' + args.End + ' ' +
            args.GenomeDir + '/ref_annot.gtf.added ' + args.Mappability + ' ' + args.GenomeDir + '/exon_probe.hg19.gene.new.bed ' + programdir + '/bin/')
    if aa == 0:
        print('Finish Read Processing! Index: ' + args.Begin + ' ~ ' + args.End)
    else:
        print('ERROR!!!!!')

if sys.argv[1] == 'FusionCandidate':
    if not args.OutDir or not args.Begin or not args.End or not args.GenomeDir:
        parser7.print_help()
        print('Please specify all required parameters!')
        exit()
    aa = os.system(
        'sh ' + programdir + '/bin/CombinePipeline_startwith_FS.sh ' + args.OutDir + ' ' + args.Begin + ' ' + args.End + ' ' + args.Prefix + ' ' + args.GenomeDir + '/ref.fa ' + args.GenomeDir + '/ref_annot.gtf.added ' + programdir + '/bin/')
    if aa == 0:
        print('Finish identifying fusion candidates')
    else:
        print('ERROR!!!!!')

if sys.argv[1] == 'Retrain':
    if not args.OutDir:
        parser8.print_help()
        print('Please specify all required parameters!')
        exit()
    aa = os.system(
        'sh ' + programdir + '/bin/CombinePipeline_Retrain.sh ' + args.OutDir + ' ' + args.Prefix + ' ' + args.Weight + ' ' + args.Epoch + ' ' + programdir + '/bin/')
    if aa == 0:
        print('Finish Retraining!')
        print('New weight file was saved at ' + args.OutDir + '/weights/RetrainWeight.hdf5')
    else:
        print('ERROR!!!!!')

if sys.argv[1] == 'ArtifactScoring':
    if not args.OutDir or not args.Begin or not args.End or not args.GenomeDir:
        parser9.print_help()
        print('Please specify all required parameters!')
        exit()
    aa = os.system(
        'sh ' + programdir + '/bin/CombinePipeline_Predict.sh ' + args.OutDir + ' ' + args.Begin + ' ' + args.End + ' ' + args.Prefix + ' ' +
        args.Weight + ' ' + args.GenomeDir + '/ref.fa ' + args.GenomeDir + '/ref_annot.gtf.added ' + programdir + '/bin/')
    if aa == 0:
        print('Finish ArtifactScoring Step!')
    else:
        print('ERROR!!!!!')

if sys.argv[1] == 'FusionReport':
    if not args.FileDir or not args.OutDir or not args.Begin or not args.End or not args.GenomeDir:
        parser10.print_help()
        print('Please specify all required parameters!')
        exit()
    numcell = 0
    for i in range(int(args.Begin), int(args.End) + 1):
        if os.path.exists(args.FileDir + str(i) + '_2.fastq') or os.path.exists(args.FileDir + str(i) + '_2.fastq.gz'):
            numcell += 1
    if args.LncRNAFilterOff:
        lncfilter = '1'
    else:
        lncfilter = '0'
    if args.NoApprovedSymbolFilterOff:
        nasfilter = '1'
    else:
        nasfilter = '0'
    aa = os.system(
        'sh ' + programdir + '/bin/CombinePipeline_startwith_ChiDist.sh ' + args.OutDir + ' ' + args.Prefix + ' ' + str(numcell) +
        ' ' + args.PvalueCutoff + ' ' + args.ArtifactScoreCutoff + ' ' + args.GenomeDir + '/ref_annot.gtf.added ' + args.GenomeDir + '/GenePos.txt ' + programdir + '/bin/ ' + lncfilter + ' ' + nasfilter)
    if aa == 0:
        if args.Prefix == '.':
            args.Prefix = ''
        print('Final Results are in ' + args.OutDir + '/FinalResult/' + args.Prefix + 'FinalOutput.abridged.txt')
        os.system(
            'cat ' + args.OutDir + '/FinalResult/' + args.Prefix + 'FinalOutput.abridged.txt > ' + args.OutDir + '/Result.abridged.txt')
        os.system(
            'cat ' + args.OutDir + '/FinalResult/' + args.Prefix + 'FinalOutput.full.txt > ' + args.OutDir + '/Result.full.txt')
        os.system('mv ' + args.OutDir + '/FinalResult/ ' + args.OutDir + '/Resulttemp/')
    else:
        print('ERROR!!!!!')

if sys.argv[1] == 'DeleteTempFiles':
    if not args.OutDir:
        parser11.print_help()
        print('Please specify all required parameters!')
        exit()
    if args.AllTempFiles:
        os.system('rm -r ' + args.OutDir + '/weights/')
        os.system('rm -r ' + args.OutDir + '/Retrain/')
        os.system('rm -r ' + args.OutDir + '/Expr/')
        os.system('rm -r ' + args.OutDir + '/ChimericOut/')
        os.system('rm -r ' + args.OutDir + '/STARMapping/')
        os.system('rm -r ' + args.OutDir + '/Resulttemp/')
        os.system('rm -r ' + args.OutDir + '/ChiDist/')
    else:
        if args.STARMappingFiles:
            os.system('rm -r ' + args.OutDir + '/STARMapping/')
        if args.AllUnimportantTempFiles:
            os.system('rm -r ' + args.OutDir + '/weights/')
            os.system('rm -r ' + args.OutDir + '/Retrain/')
            os.system('rm -r ' + args.OutDir + '/Expr/')
            os.system('rm -r ' + args.OutDir + '/ChimericOut/')
            os.system('rm -r ' + args.OutDir + '/Resulttemp/*.*')
            os.system('rm -r ' + args.OutDir + '/ChiDist/*.npy')
            os.system('rm -r ' + args.OutDir + '/ChiDist/ChiDist_middle.txt')
            os.system('rm -r ' + args.OutDir + '/ChiDist/FusionRead.txt')
            os.system('rm -r ' + args.OutDir + '/ChiDist/Homo.txt')
            os.system('rm -r ' + args.OutDir + '/ChiDist/Prob.txt')
            os.system('rm -r ' + args.OutDir + '/ChiDist/ChiDist_filtered.txt')
    print('Finish!')



