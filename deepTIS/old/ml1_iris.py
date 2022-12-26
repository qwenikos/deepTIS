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

#################  END GLOBAL  ################

#################  FUNCTIONS  #################

def createTrainingSetFunc():
    import numpy as np
    from sklearn.datasets import load_iris
    from keras.utils import to_categorical
    from sklearn.model_selection import train_test_split
    

    printFuncName()
    printn("start")
    iris =load_iris()
    type (iris)
    # sklearn.utils.Bunch
    # print (iris.DESCR)
    x=iris.data
    print(x)
    y=iris.target
    # print(y)
    #  conver to one hot encoding
    y=to_categorical(y)
    print(y.shape)
    # print (y)
    xTrain,xTest,yTrain,yTest=train_test_split(x,y,test_size=0.33,random_state=42)


    return xTrain,xTest,yTrain,yTest


def trainingFunc(xTrain,yTrain):
    printFuncName()
    from sklearn.preprocessing import MinMaxScaler
    from keras.models import Sequential
    from keras.layers import Dense 

    scaler_object=MinMaxScaler()
    scaler_object.fit(xTrain)
    xTrainScaled=scaler_object.transform(xTrain)
    
    model=Sequential()
    model.add(Dense(8,input_dim=4,activation='relu'))
    model.add(Dense(8,input_dim=4,activation='relu'))
    model.add(Dense(3,input_dim=4,activation='softmax')) ##for propabilities
    model.compile(loss='categorical_crossentropy',optimizer='adam',metrics=['accuracy'])
    print(model.summary())
    # model.fit(xTrainScaled,yTrain,epochs=150,verbose=1)

    model.fit(xTrainScaled,yTrain,epochs=150,verbose=0)
    model.save(modelName)


def predictFunc(xTest):
    from sklearn.preprocessing import MinMaxScaler
    printFuncName()
    from keras.models import load_model
    model=load_model(modelName)
    scaler_object=MinMaxScaler()
    scaler_object.fit(xTest)
    xTestScaled=scaler_object.transform(xTest)
    predictions=model.predict(xTestScaled)
    predictedClasses=np.argmax(predictions,1)
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

seqIn="ACATGAAAAAAACCCCC"
startingPos=6
length=7
extractedSeq=seqIn[startingPos-1:startingPos-1+length]
print (extractedSeq)

encoded=oneHotEncoder(seqIn)
print (encoded)


# xTrain,xTest,yTrain,yTest=createTrainingSetFunc()
# trainingFunc(xTrain,yTrain)

# predictions,predictedClasses=predictFunc(xTest)
# predictReport(predictedClasses,yTest)

# visualization
modelVisualization()



