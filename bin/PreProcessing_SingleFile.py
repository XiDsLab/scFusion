from keras.models import Sequential
from keras.layers import Embedding,Dropout,Bidirectional,Flatten,Dense,LSTM,TimeDistributed
from keras.callbacks import ModelCheckpoint,CSVLogger
from keras.layers import Conv1D, GlobalAveragePooling1D, MaxPooling1D
import numpy as np
import sys
########################################################################################################################
def ReverseComplement(str):
    return str[::-1].replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g').upper()



np.random.seed(1122)
readfile = sys.argv[1]
findgang = sys.argv[1].rfind('/')
filedir = sys.argv[1][:findgang+1]
outprefix = ''
if len(sys.argv) == 3:
    outprefix = sys.argv[2]

with open(readfile,'r') as f:
    ChimericRead_info = f.read()

ChimericRead_info = ChimericRead_info.split('\n')

ChimericRead =[]
ChimericRead_rev = []
# ChimericPoint=[]
Cont = 1
for readinfo in ChimericRead_info:
    readinfo_split = readinfo.split('\t')
    if len(readinfo_split) <= 1:
        continue
    read = readinfo_split[0]
    MergePoint = int(readinfo_split[1])
    read_new = read[0:MergePoint]+'H'+read[MergePoint:]
    # read_new = read
    Point = np.zeros(60)
    Point[MergePoint-1] = 1
    Point[MergePoint] = 1
    if 'N' not in read_new:
        ChimericRead.append(read_new)
        ChimericRead_rev.append(ReverseComplement(read_new))
    Cont = Cont + 1
########################################################################################################################

Data1 = np.ndarray(shape=(len(ChimericRead),61,1),dtype=float)
Data2 = np.ndarray(shape=(len(ChimericRead_rev),61,1),dtype=float)
for index in range(len(ChimericRead)):
    Data1[index,:,0] = np.array([int(c) for c in ChimericRead[index].upper().replace('A','0').replace('T','1').replace('C','2').replace('G','3').replace('H','4')])
    Data2[index,:,0] = np.array([int(c) for c in ChimericRead_rev[index].upper().replace('A','0').replace('T','1').replace('C','2').replace('G','3').replace('H','4')])

    
########################################################################################################################

LIST1 = list(range(len(ChimericRead)))
LIST2 = list(range(len(ChimericRead_rev)))

DataNum = len(ChimericRead)
Good_for_Tra = Data1[LIST1[0:DataNum],:,:]
Good_for_Tra_rev = Data2[LIST2[0:DataNum],:,:]

np.save(filedir + outprefix + 'Reads.npy',Good_for_Tra)
np.save(filedir + outprefix + 'Reads_rev.npy',Good_for_Tra_rev)