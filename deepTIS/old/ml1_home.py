# ml1.py
from lib.ml1_lib import *
from lib.general_lib import *
import os
import numpy as np

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'   ## to hide some tensorFlow warnings

#################  GLOBAL  ####################

modelDir="models/"
imgDir="img/"
os.makedirs(modelDir, exist_ok=True)
os.makedirs(imgDir, exist_ok=True)

modelName=modelDir+"TIS_model.nic"

workingDrir="/home/nikos/MEGAsync/UOA/Bioinformatics/"
# workingDrir="/mnt/raid1/perdikopn/local/"
positiveSetFileName=workingDrir+"DeepTISFramework/trainingData/TIS/positive/sequences/tisSeq_valid.fa"
negativeSetFileName=workingDrir+"DeepTISFramework/trainingData/TIS/negative/sequences/tisSeq_valid.fa"

ensembleTestDatasetFileName=workingDrir+"DeepTISFramework/testData/sorf.org/allSorfsWithTIS_0_100_ENSEMBLE_only.fa"
intronicTestDatasetFileName=workingDrir+"DeepTISFramework/testData/sorf.org/allSorfsWithTIS_0_100_intronic_only.fa"
pos_sim_10k=workingDrir+"DeepTISFramework/testData/simmulated_data/positive_simSorfs_Flank-100_10000.fa"
neg_sim_10k=workingDrir+"DeepTISFramework/testData/simmulated_data/negative_simSorfs_Flank-100_10000.fa"

#################  END GLOBAL  ################

#################  FUNCTIONS  #################


def readTrainingData():
    printFuncName()

    import numpy as np

    size=8000 ##8000
    ATGPos=100

    x=[]
    y=[]
    
    cnt=0
    aPositiveFile=open(positiveSetFileName,"r")
    for aLine in aPositiveFile:
        if aLine[0]==">":
            continue
        trainSeq=aLine[ATGPos-7:ATGPos]               #part befote ATG
        trainSeq+=aLine[ATGPos+3:ATGPos+10]            #part after ATG
        x+=[trainSeq]
        y+=[1]  #o for positive
        cnt+=1
        if (cnt==size):
            break

    cnt=0
    aNegativeFile=open(negativeSetFileName,"r")
    for aLine in aNegativeFile:
        if aLine[0]==">":
            continue
        trainSeq=aLine[ATGPos-7:ATGPos]               #part befote ATG
        trainSeq+=aLine[ATGPos+3:ATGPos+10]            #part after ATG
        x+=[trainSeq]
        y+=[0] #1 for negative
        cnt+=1
        if (cnt==size):
            break

    return x,y

def readTestData():
    printFuncName()
    import numpy as np

    size=8000 ##8000
    # ATGPos=100
    ATGPos=8

    x=[]
    y=[]
    
    cnt=0
    aPositiveFile=open(ensembleTestDatasetFileName,"r")
    ATGPos=8
    # aPositiveFile=open(intronicTestDatasetFileName,"r")
    # aPositiveFile=open(neg_sim_10k,"r")
    # ATGPos=100
    intronicTestDatasetFileName
    for aLine in aPositiveFile:
        if aLine[0]==">":
            continue
        trainSeq=aLine[ATGPos-7:ATGPos]               #part befote ATG
        trainSeq+=aLine[ATGPos+3:ATGPos+10]            #part after ATG
        x+=[trainSeq]
        y+=[0]  #0 for positive
        cnt+=1
        if (cnt==size):
            break
    return x,y

def createNetInput(x,y):
    printFuncName()
    import numpy as np
    netX=[]
    netY=[]

    for aSeq,y in zip(x,y):
        # print(aSeq,y)
        # encodedSeq=oneHotEncoder1D(aSeq)
        encodedSeq=oneHotEncoder(aSeq)
        netX+=[encodedSeq]
        netY+=[y]

    return np.array(netX),np.array(netY)

def splitTrainTestData(x,y):
    printFuncName()
    from sklearn.model_selection import train_test_split
    xTrain,xTest,yTrain,yTest=train_test_split(x,y,test_size=0.33,random_state=42)
    return xTrain,xTest,yTrain,yTest


def trainingFunc(xTrain,yTrain):
    printFuncName()
    # from keras.layers.convolutional import MaxPooling1D,Conv1D
    from keras.layers import Input,Conv1D,Flatten,Dense,MaxPooling1D,Dropout,LSTM,Activation
    from keras.constraints import maxnorm
    # from sklearn.preprocessing import MinMaxScaler
    from keras.models import Sequential  
    # from keras.layers.convolutional import Conv1D ,Convolution1D 
    # from keras.layers import LSTM,Dense, Flatten

    
    print (">> xTrain.shape\t",xTrain.shape)
    print (">> yTrain.shape\t",yTrain.shape)

    modelNo=2
    model=Sequential()
    if modelNo==1:
        model.add(Input(shape=(9,4)))  ##input shape
        model.add(Conv1D(128, 3, activation='relu'))
        model.add(Flatten())
        model.add(Dense(32, activation='relu'))
        model.add(Dense(1, activation='softmax')) ##output

        model.compile(optimizer='adam',loss='categorical_crossentropy',metrics=['accuracy'])
        model.summary()
        # model.fit(xTrain,yTrain, epochs=5, steps_per_epoch=100)
        model.fit(xTrain,yTrain, epochs=5, )

    if modelNo==2:
        # model.add(Conv1D(nb_filter=128,    ##number of conv kernel to use ,dim of output)
        #                     filter_length=3,      ##The extension (spatial or temporal) of each filter.
        #                     input_dim=8,          ##Number of channels/dimensions in the input. Either this argument or the keyword argument input_shape must be provided when using this layer as the first layer in a model.
        #                     input_length=203,     ##Length of input sequences, when it is constant. This argument is required if you are going to connect Flatten then Dense layers upstream (without it, the shape of the dense outputs cannot be computed).
        #                     border_mode='valid',  ##'valid', 'same' or 'full' ('full' requires the Theano backend).
        #                     W_constraint = maxnorm(3),  ## instance of the constraints module (eg. maxnorm, nonneg), applied to the main weights matrix.
        #                     activation='relu',
        #                     subsample_length=1))  ## factor by which to subsample output.

        model.add(Input(shape=(14,4)))  ##input shape
        model.add(Conv1D(filters=128,kernel_size=4, activation='relu')) #filters kernel_size
        model.add(MaxPooling1D(pool_size=3))   ##Integer, size of the max pooling window.
        model.add(Dropout(rate=0.21370950078747658))
        model.add(LSTM(units=256,return_sequences=True))
        model.add(Dropout(rate=0.7238091317104384))
        model.add(Flatten())
        model.add(Dense(1))
        model.add(Activation('sigmoid'))

        print(model.summary())

        # model.compile(loss='binary_crossentropy',
        #           optimizer='nadam',
        #           metrics=['accuracy'])
        # model.compile(optimizer='adam',loss='categorical_crossentropy',metrics=['accuracy'])
        model.compile(optimizer='adam',loss='binary_crossentropy',metrics=['accuracy'])
        model.summary()
        model.fit(xTrain,yTrain,epochs=150,verbose=1)
    model.save(modelName)


def predictFunc(xTest):
    printFuncName()
    from keras.models import load_model
    model=load_model(modelName)

    predictions=model.predict(xTest)
    predictedClasses=np.copy(predictions)
    predictedClasses[predictedClasses>=0.5]=1
    predictedClasses[predictedClasses<0.5]=0
    # predictedClasses=np.argmax(predictions,1)
    return predictions,predictedClasses


def predictReport(predictedClasses,yTest):
    printFuncName()
    from sklearn.metrics import confusion_matrix,classification_report, accuracy_score

    print (predictedClasses)
    testClasses=np.argmax(yTest,1)
    print(testClasses)
    confMat=confusion_matrix(predictedClasses,testClasses)
    print (confMat)
    print (classification_report(testClasses,predictedClasses ))


def modelVisualization():
    import visualkeras
    from keras.models import load_model
    modelImage=imgDir+"model.png"
    model=load_model(modelName)
    # visualkeras.layered_view(model, legend=True,to_file =modelImage ).show()
    visualkeras.layered_view(model,to_file =modelImage, legend=True).show

#################  END FUNCTIONS  #############

################# MAIN PROGRAM  ###############
###############################################




printn("start")

x,y=readTrainingData()

netX,netY=createNetInput(x,y)
# print (netX.shape)

xTrain,xTest,yTrain,yTest=splitTrainTestData(netX,netY)
## have to save the sets

training=True
if training:
    trainingFunc(xTrain,yTrain)
    predictions,predictedClasses=predictFunc(xTest)
    # print (predictedClasses)
    # predictReport(predictedClasses,yTest)
    exit()

xTest,yTest=readTestData()
xNet,yNet=createNetInput(xTest,yTest)
predictions,predictedClasses=predictFunc(xNet)
# print (predictedClasses)
print((predictedClasses == 0.).sum())
print((predictedClasses == 1.).sum())
# predictReport(predictedClasses,yTest)
# print (predictions)

# trainingFunc(xTrain,yTrain)

# predictions,predictedClasses=predictFunc(xTest)
# predictReport(predictedClasses,yTest)

# visualization
# modelVisualization()




