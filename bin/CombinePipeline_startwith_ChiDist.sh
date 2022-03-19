#!/bin/bash

FilePath=$1
prefix=$2
numcell=$3
pv=$4
fakeprob=$5
gtf=$6
genepos=$7
codedir=$8
lncfilter=$9
nasfilter=${10}

if [ "${prefix}" = "." ]
then
	prefix=""
fi

mkdir -p ${FilePath}/FinalResult/
mkdir -p ${FilePath}/FinalResult/temp/
Rscript ${codedir}/ChiDist_2.2.1_Combine.R ${numcell} ${FilePath}/ChiDist/${prefix}ChiDist.txt ${FilePath}/ChiDist/${prefix}ChiDist_filtered.txt ${FilePath}/FinalResult/temp/${prefix}Allresult.txt ${FilePath}/FinalResult/temp/${prefix}pars.RData
python ${codedir}/TidyupFusionFinalResult.py ${FilePath}/FinalResult/temp/${prefix}Allresult.txt ${FilePath}/ChimericOut/ > ${FilePath}/FinalResult/temp/${prefix}Allresult_filtered.txt
python ${codedir}/Results_Filtered2Final.py ${FilePath}/FinalResult/temp/${prefix}Allresult_filtered.txt ${pv} ${fakeprob} ${numcell} > ${FilePath}/FinalResult/${prefix}Final.txt
python ${codedir}/TidyupFusionFinalResult_FindSupCell.py ${FilePath}/FinalResult/${prefix}Final.txt ${FilePath}/ChimericOut/ > ${FilePath}/FinalResult/${prefix}Final_Cells.txt
python ${codedir}/ResultLastFiltered.py ${FilePath}/FinalResult/${prefix}Final_Cells.txt ${genepos} ${FilePath}/FinalResult/${prefix}Final_Cells_filtered.txt ${gtf} ${genepos} ${lncfilter} ${nasfilter}
python ${codedir}/ResultFinalOutput.py ${FilePath}/FinalResult/${prefix}Final_Cells_filtered.txt ${gtf} ${FilePath}/ChimericOut/ ${FilePath}/FinalResult/${prefix}FinalOutput
