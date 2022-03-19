from keras.models import Sequential
from keras.layers import Embedding,Dropout,Bidirectional,Flatten,Dense,LSTM,TimeDistributed
from keras.callbacks import ModelCheckpoint,CSVLogger
from keras.layers import Conv1D, GlobalAveragePooling1D, MaxPooling1D
import numpy as np
import sys
########################################################################################################################
np.random.seed(1122)
ChimericFile = sys.argv[1]
FakeChimericFile = sys.argv[2]
Outdir = sys.argv[3]
with open(ChimericFile,'r') as f:
    ChimericRead_info = f.read()

ChimericRead_info = ChimericRead_info.split('\n')

ChimericRead =[]
# ChimericPoint=[]
Cont = 1
for readinfo in ChimericRead_info:
    readinfo_split = readinfo.split('\t')
    read = readinfo_split[0]
    try:
        MergePoint = int(readinfo_split[1])
    except:
        sys.stderr.write(readinfo)
    read_new = read[0:MergePoint]+'H'+read[MergePoint:]
    # read_new = read
    Point = np.zeros(60)
    Point[MergePoint-1] = 1
    Point[MergePoint] = 1
    if 'N' not in read_new:
        ChimericRead.append(read_new)
        read_new_inv = read_new.replace('A','5').replace('T','6').replace('C','7').replace('G','8')
        read_new_inv = read_new_inv[::-1]
        read_new_inv = read_new_inv.replace('5','T').replace('6','A').replace('7','G').replace('8','C')
        ChimericRead.append(read_new_inv)
        # ChimericPoint.append(Point)
    # print(Cont)
    Cont = Cont + 1
########################################################################################################################
with open(FakeChimericFile,'r') as f:
    ManmadeRead_info = f.read()

ManmadeRead_info = ManmadeRead_info.split('\n')

ManmadeRead =[]
# ManmadePoint=[]
Cont = 1
for readinfo in ManmadeRead_info:
    readinfo_split = readinfo.split('\t')
    read = readinfo_split[0]
    try:
        MergePoint = int(readinfo_split[1])
    except:
        sys.stderr.write(readinfo)
    read_new = read[0:MergePoint]+'H'+read[MergePoint:]
    # read_new = read
    Point = np.zeros(61)
    Point[MergePoint-1] = 1
    Point[MergePoint] = 1

    if 'N' not in read_new:
        ManmadeRead.append(read_new)
        read_new_inv = read_new.replace('A','5').replace('T','6').replace('C','7').replace('G','8')
        read_new_inv = read_new_inv[::-1]
        read_new_inv = read_new_inv.replace('5','T').replace('6','A').replace('7','G').replace('8','C') 
        ManmadeRead.append(read_new_inv)
        # ManmadePoint.append(Point)

    # print(Cont)
    Cont = Cont + 1
########################################################################################################################

Data1 = np.ndarray(shape=(len(ChimericRead),61,1),dtype=float)
for index in range(len(ChimericRead)):
    Data1[index,:,0] = np.array([int(c) for c in ChimericRead[index].replace('A','0').replace('T','1').replace('C','2').replace('G','3').replace('H','4')])
    # Data1[index,:,1] = ManmadePoint[index]

Data2 = np.ndarray(shape=(len(ManmadeRead),61,1),dtype=float)
for index in range(len(ManmadeRead)):
    Data2[index,:,0] = np.array([int(c) for c in ManmadeRead[index].replace('A','0').replace('T','1').replace('C','2').replace('G','3').replace('H','4')])
    # Data2[index,:,1] = ManmadePoint[index]

########################################################################################################################

DataNum = min(Data1.shape[0],Data2.shape[0])
Data1 = Data1[0:DataNum,:,:]
Data2 = Data2[0:DataNum,:,:]

TraNum = int(DataNum*0.7)

Good_for_Tra = Data1[0:TraNum,:,:]
Simu_for_Tra = Data2[0:TraNum,:,:]

Good_for_Tst = Data1[TraNum:,:,:]
Simu_for_Tst = Data2[TraNum:,:,:]


LIST = list(range(Good_for_Tra.shape[0]))
np.random.shuffle(LIST)
Good_for_Tra = Good_for_Tra[LIST,:,:]
Simu_for_Tra = Simu_for_Tra[LIST,:,:]


LIST = list(range(Good_for_Tst.shape[0]))
np.random.shuffle(LIST)
Good_for_Tst = Good_for_Tst[LIST,:,:]
Simu_for_Tst = Simu_for_Tst[LIST,:,:]
    

np.save(Outdir + '/Good_for_Tra.npy',Good_for_Tra)
np.save(Outdir + '/Simu_for_Tra.npy',Simu_for_Tra)

np.save(Outdir + '/Good_for_Tst.npy',Good_for_Tst)
np.save(Outdir + '/Simu_for_Tst.npy',Simu_for_Tst)

########################################################################################################################





