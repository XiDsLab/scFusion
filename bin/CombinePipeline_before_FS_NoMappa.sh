#!/bin/bash

FilePath=$1
mystart=$2
myend=$3
gtffile=$4
exonfile=$5
codedir=$6
mkdir -p ${FilePath}/STARMapping
mkdir -p ${FilePath}/ChimericOut
mkdir -p ${FilePath}/Expr/


for ((i=${mystart};i<=${myend};i++))
do
	file=`ls ${FilePath}/STARMapping/${i}/*Chimeric.out.sam`
	if [[ -n ${file} ]]; then
		python ${codedir}/RmLowMappibility_ChimericRead_NoFilter.py ${file} ${FilePath}/ChimericOut/${i}.sam
	fi
done

sh ${codedir}/Annotate.sh ${FilePath}/ChimericOut/ ${mystart} ${myend} ${gtffile} ${codedir}
sh ${codedir}/FindFusionSupport.sh ${FilePath}/ChimericOut/ ${mystart} ${myend} ${codedir}
sh ${codedir}/CalcRPKM.sh ${exonfile} ${FilePath}/STARMapping ${FilePath}/Expr/ ${mystart} ${myend}