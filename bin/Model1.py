# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 21:38:29 2019

@author: BioMed-X
"""


from keras.models import Sequential,Model
from keras.layers import Embedding,Dropout,Bidirectional,Flatten,Dense,LSTM,TimeDistributed, Activation,Input,merge,concatenate
from keras.callbacks import ModelCheckpoint,CSVLogger
from keras.layers import Conv1D, GlobalAveragePooling1D, MaxPooling1D
from tensorflow.keras.optimizers import Adam
import numpy as np
from tensorflow.keras.utils import to_categorical
import os

os.environ["CUDA_VISIBLE_DEVICES"] = "0"

def Cla_LSTM():
# 搭模型
    INPUT1 = Input(shape=(61,))
    # INPUT2 = Input(shape=(60,1))
    INPUT1_enco = Embedding(5,5,input_length=61)(INPUT1)

    # MERGE = merge((INPUT1_enco, INPUT2),mode='concat',concat_axis=-1)
    # MERGE = concatenate([INPUT1,INPUT2])
    LSTM1 = Bidirectional(LSTM(32,return_sequences=True),merge_mode='concat')(INPUT1_enco)
    DROP1 = Dropout(0.5)(LSTM1)
    
    LSTM2 = Bidirectional(LSTM(64,return_sequences=True),merge_mode='concat')(DROP1)
    DROP2 = Dropout(0.5)(LSTM2)
    
    LSTM3 = Bidirectional(LSTM(128,return_sequences=True),merge_mode='concat')(DROP2)
    DROP3 = Dropout(0.5)(LSTM3)
    
    '''
    LSTM4 = Bidirectional(LSTM(64,return_sequences=True),merge_mode='concat')(DROP3)
    DROP4 = Dropout(0.5)(LSTM4)    
    
    LSTM5 = Bidirectional(LSTM(128,return_sequences=True),merge_mode='concat')(DROP4)
    DROP5 = Dropout(0.5)(LSTM5) 

    LSTM6 = Bidirectional(LSTM(128,return_sequences=True),merge_mode='concat')(DROP5)
    DROP6 = Dropout(0.5)(LSTM6) 

    LSTM7 = Bidirectional(LSTM(36,return_sequences=True),merge_mode='concat')(DROP6)
    DROP7 = Dropout(0.5)(LSTM7) 
    '''
    LSTM8 = Bidirectional(LSTM(256, return_sequences=False), merge_mode='concat')(DROP3)
    DENSE1 = Dense(256)(LSTM8)
    DENSE2 = Dense(2)(DENSE1)
    
    ACT1 = Activation('softmax')(DENSE2)
    # model = Model(inputs=[INPUT1,INPUT2],outputs= ACT1)
    model = Model(inputs=INPUT1,outputs= ACT1)
    return model

if __name__ == '__main__':
    np.random.seed(1122)
    Good_for_Tra = np.load('Good_for_Tra.npy')
    Simu_for_Tra = np.load('Simu_for_Tra.npy')
    Good_for_Tst = np.load('Good_for_Tst.npy')
    Simu_for_Tst = np.load('Simu_for_Tst.npy')
    Tra_x = np.squeeze(np.concatenate((Good_for_Tra,Simu_for_Tra),axis=0))
    Tra_y = np.concatenate( (np.zeros((Good_for_Tra.shape[0],1)),np.ones((Simu_for_Tra.shape[0],1))), axis=0)
    Tst_x = np.squeeze(np.concatenate((Good_for_Tst,Simu_for_Tst),axis=0))
    Tst_y = np.concatenate((np.zeros((Good_for_Tst.shape[0], 1)), np.ones((Simu_for_Tst.shape[0], 1))), axis=0)
    Tra_y = to_categorical(Tra_y)
    Tst_y = to_categorical(Tst_y)

    LIST = list(range(Tra_x.shape[0]))
    np.random.shuffle(LIST)

    Tra_x = Tra_x[LIST,:]
    Tra_y = Tra_y[LIST,:]
    # Tra_x_input1 = Tra_x[..., 0]
    # Tra_x_input2 = Tra_x[..., 1][...,np.newaxis]

    LIST = list(range(Tst_x.shape[0]))
    np.random.shuffle(LIST)
    Tst_x = Tst_x[LIST,:]
    Tst_y = Tst_y[LIST,:]

    model = Cla_LSTM()
    model.load_weights('weight-010.hdf5')

    ADAM = Adam(learning_rate=0.001)
    model_checkpoint = ModelCheckpoint(filepath='weight-{epoch:03d}.hdf5', verbose=1, monitor='val_loss', save_best_only=True)

    model.compile(loss='binary_crossentropy', optimizer=ADAM, metrics=['accuracy'])
    csv_loger=CSVLogger('log.csv',append=True,separator=';')

    # 训练模型
    batch_size = 2500
    epochs = 200
    
    # model.fit(x=[Tra_x_input1,Tra_x_input2], y=Tra_y,batch_size=batch_size,epochs=epochs, verbose=1 ,callbacks=[model_checkpoint,csv_loger], validation_split=0.25, shuffle=True)
    model.fit(x=Tra_x, y=Tra_y,batch_size=batch_size,epochs=epochs, initial_epoch= 11,validation_data=(Tst_x, Tst_y),verbose=1 ,callbacks=[model_checkpoint,csv_loger])


# 

'''
model.fit(x_train, y_train, epochs=20, batch_size=128)
score = model.evaluate(x_test, y_test, batch_size=128)


model.add(Conv1D(filters=2, kernel_size=1 , activation='relu', name='conv1_1'))

model.add(Conv1D(filters=64, kernel_size=5 , activation='relu', name='conv1_1'))
model.add(MaxPooling1D(2))

model.add(Conv1D(filters=64, kernel_size=5 , activation='relu', name='conv1_1'))

model.add(Conv1D(filters=2, kernel_size=2 , activation='relu', name='conv1_1'))

model.add(Bidirectional(LSTM(20,return_sequences=True),merge_mode='concat'))

model.summary()

model.add(Activation('softmax'));
'''

    