#!/bin/bash

FilePath=$1
mystart=$2
myend=$3
prefix=$4
hg19file=$5
gtf=$6
codedir=$7

if [ "${prefix}" = "." ]
then
	prefix=""
fi

mkdir -p ${FilePath}/ChiDist/
python ${codedir}/FusionScore.py ${FilePath}/ChimericOut/ ${mystart} ${myend} ${FilePath}/Expr/ > ${FilePath}/ChimericOut/${prefix}FusionScore.txt
python ${codedir}/FindHomoPattern_RAM.py ${FilePath}/ChimericOut/${prefix}FusionScore.txt ${hg19file} ${gtf} > ${FilePath}/ChiDist/${prefix}Homo.txt
python ${codedir}/FindChiDist.py ${FilePath}/ChimericOut/ ${mystart} ${myend} ${FilePath}/Expr/ ${FilePath}/ChiDist/${prefix}Homo.txt ${prefix} > ${FilePath}/ChiDist/${prefix}ChiDist_middle.txt
python ${codedir}/PreProcessing_SingleFile.py ${FilePath}/ChiDist/${prefix}FusionRead.txt ${prefix}
