# -*- coding: utf-8 -*-
"""
@author: Huang Wenjian
"""

from keras.models import Sequential
from keras.layers import Embedding,Dropout,Bidirectional,Flatten,Dense,LSTM,TimeDistributed, Activation
from keras.callbacks import ModelCheckpoint,CSVLogger
from keras.layers import Conv1D, GlobalAveragePooling1D, MaxPooling1D
from tensorflow.keras.optimizers import Adam
import numpy as np
from tensorflow.keras.utils import to_categorical
from Model1 import Cla_LSTM
import os
import sys

# os.environ["CUDA_VISIBLE_DEVICES"] = "0"
np.random.seed(1122)
outfile = open(sys.argv[1], 'w')
weightfile = sys.argv[2]
prefix = ''
if len(sys.argv) == 4:
    prefix = sys.argv[3]

findgang = sys.argv[1].rfind('/')
filedir = sys.argv[1][:findgang+1]
findgang = sys.argv[0].rfind('/')
codedir = sys.argv[0][:findgang+1]

Good_for_Tra = np.load(filedir + prefix + 'Reads.npy')
Good_for_Tra_rev = np.load(filedir + prefix + 'Reads_rev.npy')
Tst_x = np.squeeze(Good_for_Tra)
Tst_x_rev = np.squeeze(Good_for_Tra_rev)

model = Cla_LSTM()
model.load_weights(weightfile)

batch_size = 500 

Prob = model.predict(Tst_x,batch_size)
Prob_rev = model.predict(Tst_x_rev,batch_size)
AveProb = (Prob[:,0] + Prob_rev[:,0]) / 2
for i in range(len(AveProb)):
    outfile.write(str(Prob[i,0]) + '\t' + str(Prob_rev[i,0]) + '\t' + str(AveProb[i]) + '\n')
outfile.close()
