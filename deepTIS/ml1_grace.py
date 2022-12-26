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
workingDrir="/mnt/raid1/perdikopn/local/"
positiveSetFileName=workingDrir+"DeepTISFramework/trainingData/TIS/positive/sequences/tisSeq_valid.fa"
negativeSetFileName=workingDrir+"DeepTISFramework/trainingData/TIS/negative/sequences/tisSeq_valid.fa"

ensembleTestDatasetFileName=workingDrir+"DeepTISFramework/testData/sorf.org/allSorfsWithTIS_0_100_ENSEMBLE_only.fa"
intronicTestDatasetFileName=workingDrir+"DeepTISFramework/testData/sorf.org/allSorfsWithTIS_0_100_intronic_only.fa"
pos_sim_10k=workingDrir+"DeepTISFramework/testData/simmulated_data/positive_simSorfs_Flank-100_10000.fa"
neg_sim_10k=workingDrir+"DeepTISFramework/testData/simmulated_data/negative_simSorfs_Flank-100_10000.fa"

upSeqSize=7
downSeqSize=30 ##not count ATG
seqSizeWithStartCodon=upSeqSize+3+downSeqSize
seqSize=upSeqSize+downSeqSize


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
        trainSeq=aLine[ATGPos-upSeqSize:ATGPos]               #part befote ATG
        trainSeq+=aLine[ATGPos+3:ATGPos+3+downSeqSize]            #part after ATG
        
        if len(trainSeq)<seqSize:
            print("Warning:length Problem in",trainSeq)
            continue
        if (dnaCheck(trainSeq)==False):
            print("Warning:Nucleotide Problem in",trainSeq)
            continue
        x+=[trainSeq]    
        y+=[1]  #0 for positive
        cnt+=1
        if (cnt==size):
            break
    cnt=0

    aNegativeFile=open(negativeSetFileName,"r")
    for aLine in aNegativeFile:
        if aLine[0]==">":
            continue
        trainSeq=aLine[ATGPos-upSeqSize:ATGPos]               #part befote ATG
        trainSeq+=aLine[ATGPos+3:ATGPos+3+downSeqSize]            #part after ATG
        
        if len(trainSeq)<seqSize:
            print("Warning:length Problem in",trainSeq)
            continue
        if (dnaCheck(trainSeq)==False):
            print("Warning:Nucleotide Problem in",trainSeq)
            continue
        x+=[trainSeq]
        y+=[0] #1 for negative
        cnt+=1
        if (cnt==size):
            break
    return x,y

def readTestData(testDataFileName,ATGPos):
    printFuncName()
    import numpy as np

    size=8000 ##8000
    
    x=[]
    y=[]
    cnt=0
    testDataFile=open(testDataFileName,"r")
    for aLine in testDataFile:
        if aLine[0]==">":
            continue
        trainSeq=aLine[ATGPos-upSeqSize:ATGPos]                   #part befote ATG
        trainSeq+=aLine[ATGPos+3:ATGPos+3+downSeqSize]            #part after ATG
        if len(trainSeq)<seqSize:
            print("Warning:length Problem in",trainSeq)
            continue
        if (dnaCheck(trainSeq)==False):
            print("Warning:Nucleotide Problem in",trainSeq)
            continue
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
        # print (aSeq)
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
    from keras.layers import Input,Conv1D,Flatten,Dense,MaxPooling1D,Dropout,LSTM,Activation
    from keras.constraints import maxnorm
    from keras.models import Sequential  
    
    inputShape=xTrain.shape[1:]
    print (">> xTrain.shape\t",xTrain.shape)
    print (">> yTrain.shape\t",yTrain.shape)

    modelNo=4
    model=Sequential()
    if modelNo==1:
        model.add(Input(shape=inputShape))  ##input shape
        model.add(Conv1D(128, 3, activation='relu'))
        model.add(Flatten())
        model.add(Dense(32, activation='relu'))
        model.add(Dense(1, activation='softmax')) ##output

        model.compile(optimizer='adam',loss='categorical_crossentropy',metrics=['accuracy'])
        model.summary()
        # model.fit(xTrain,yTrain, epochs=5, steps_per_epoch=100)
        model.fit(xTrain,yTrain, epochs=5, )

    if modelNo==2:
        model.add(Input(shape=inputShape))  ##input shape
        model.add(Conv1D(filters=128,kernel_size=4, activation='relu')) #filters kernel_size
        model.add(MaxPooling1D(pool_size=3))   ##Integer, size of the max pooling window.
        model.add(Dropout(rate=0.21370950078747658))
        model.add(LSTM(units=256,return_sequences=True))
        model.add(Dropout(rate=0.7238091317104384))
        model.add(Flatten())
        model.add(Dense(1))
        model.add(Activation('sigmoid'))
        # print(model.summary())
        # model.compile(optimizer='adam',loss='categorical_crossentropy',metrics=['accuracy'])
        model.compile(optimizer='adam',loss='binary_crossentropy',metrics=['accuracy','mse']) ## loss='categorical_crossentropy' in case of multiclass
        model.summary()
        model.fit(xTrain,yTrain,epochs=150,verbose=1)

    if modelNo==3: #https://www.nature.com/articles/s41598-021-89850-9
        model.add(Input(shape=inputShape))  ##input shape
        model.add(Conv1D(filters=128,kernel_size=3, activation='relu')) #filters kernel_size
        model.add(Conv1D(filters=64,kernel_size=3, activation='relu')) #filters kernel_size
        model.add(Flatten())
        model.add(Dense(1))
        model.add(Activation('sigmoid'))
        print(model.summary())
        model.compile(optimizer='adam',loss='binary_crossentropy',metrics=['accuracy','mse']) ## loss='categorical_crossentropy' in case of multiclass
        model.summary()
        model.fit(xTrain,yTrain,epochs=50,verbose=1)

    if modelNo==4: #nikos1
        model.add(Input(shape=inputShape))  ##input shape
        model.add(Conv1D(filters=256,kernel_size=4, activation='relu')) #filters kernel_size
        model.add(Conv1D(filters=64,kernel_size=3, activation='relu')) #filters kernel_size
        model.add(Conv1D(filters=16,kernel_size=3, activation='relu')) #filters kernel_size
        model.add(Flatten())
        model.add(Dense(1))
        model.add(Activation('sigmoid'))
        print(model.summary())
        model.compile(optimizer='adam',loss='binary_crossentropy',metrics=['accuracy','mse']) ## loss='categorical_crossentropy' in case of multiclass
        model.summary()
        model.fit(xTrain,yTrain,epochs=3,verbose=2)
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

def main():
    

    training=True
    testing=True

    if training:
        printn("start Training")
        x,y=readTrainingData()
        netX,netY=createNetInput(x,y)
        xTrain,xTest,yTrain,yTest=splitTrainTestData(netX,netY)
        ## have to save the sets in numpy array
        trainingFunc(xTrain,yTrain)
        predictions,predictedClasses=predictFunc(xTest)
        # print (predictedClasses)
        # predictReport(predictedClasses,yTest)
    
    if (testing):

        printh("ensembleTestDatasetFileName")
        xTest,yTest=readTestData(ensembleTestDatasetFileName,8)
        xNet,yNet=createNetInput(xTest,yTest)
        predictions,predictedClasses=predictFunc(xNet)
        # print (predictedClasses)
        print((predictedClasses == 0.).sum())
        print((predictedClasses == 1.).sum())

        printh("intronicTestDatasetFileName")
        xTest,yTest=readTestData(intronicTestDatasetFileName,8)

        xNet,yNet=createNetInput(xTest,yTest)
        predictions,predictedClasses=predictFunc(xNet)
        # print (predictedClasses)
        print((predictedClasses == 0.).sum())
        print((predictedClasses == 1.).sum())

        printh("pos_sim_10k")
        xTest,yTest=readTestData(pos_sim_10k,100)
        xNet,yNet=createNetInput(xTest,yTest)
        predictions,predictedClasses=predictFunc(xNet)
        # print (predictedClasses)
        print((predictedClasses == 0.).sum())
        print((predictedClasses == 1.).sum())

        printh("neg_sim_10k")
        xTest,yTest=readTestData(neg_sim_10k,100)
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



if __name__ == "__main__":
    main()