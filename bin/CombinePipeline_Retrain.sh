
DatasetName=$1
prefix=$2
initweight=$3
epochnum=$4
codedir=$5
mkdir -p ${DatasetName}/Retrain
mkdir -p ${DatasetName}/weights/

if [ "${prefix}" = "." ]
then
	prefix=""
fi

python ${codedir}/ExtractChimericRead4Retrain.py ${DatasetName}/ChiDist/${prefix}ChiDist_middle.txt > ${DatasetName}/Retrain/${prefix}ChimericRead.txt
python ${codedir}/ExtractSimulatedChimericRead4Retrain.py ${DatasetName}/Retrain/${prefix}ChimericRead.txt ${DatasetName}/STARMapping/ > ${DatasetName}/Retrain/${prefix}SimuRead.txt
python ${codedir}/Data_preprocess_MyRetrain.py ${DatasetName}/Retrain/${prefix}ChimericRead.txt ${DatasetName}/Retrain/${prefix}SimuRead.txt ${DatasetName}/Retrain/
python ${codedir}/Model1_Retrain.py ${DatasetName}/Retrain/ ${initweight} ${DatasetName}/weights/ ${epochnum}