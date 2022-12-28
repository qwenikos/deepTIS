def printn(str):
    print("\n-------  ",str,"  ------\n")

def printh(str):
    newStr=" "*int((70-len(str))/2)+str+" "*int((70-len(str))/2)
    strLen=len(newStr)
    print("================================================================================")
    print("=====                                                                     ======")
    print("====="+newStr+"="*(80-strLen-5))
    print("=====                                                                     ======")
    print("================================================================================")

def printFuncName(commentFlag):
    ### return name of function where prinfFuncName called
    ### not show commpent if only commentsFlag=0

    if commentFlag==0:
        return
    else:
        import sys
        import inspect
        functionNameAsString = inspect.stack()[1][3]
        print ("\n_________________________________________________________________")
        print ("    >>>>>>> Function Name = ",functionNameAsString,"     <<<<< " ) 
        print ("_________________________________________________________________\n")
        print(commentFlag)

def oneHotEncoder(seqIn):
    import numpy as np
    encSeq=[]
    seqIn=seqIn.upper()
    # print (seqIn)
    seqToOneHotDict={"A":[1.0, 0.0, 0.0, 0.0],"C":[0.0, 1.0, 0.0, 0.0],"T":[0.0, 0.0, 1.0, 0.0],"G":[0.0, 0.0, 0.0, 1.0]}
    for aBase in seqIn:
        encSeq+=[seqToOneHotDict[aBase]]  
    encSeq=np.array(encSeq)
    return encSeq

def oneHotEncoder1D(seqIn):
    import numpy as np
    encSeq=[]
    seqIn=seqIn.upper()
    # print (seqIn)
    seqToOneHotDict={"A":[1.0, 0.0, 0.0, 0.0],"C":[0.0, 1.0, 0.0, 0.0],"T":[0.0, 0.0, 1.0, 0.0],"G":[0.0, 0.0, 0.0, 1.0]}
    for aBase in seqIn:
        encSeq+=seqToOneHotDict[aBase]
    encSeq=np.array(encSeq)
    return encSeq


def oneHotEncoder1D(seqIn):
    import numpy as np
    encSeq=[]
    seqIn=seqIn.upper()
    # print (seqIn)
    seqToOneHotDict={"A":[1.0, 0.0, 0.0, 0.0],"C":[0.0, 1.0, 0.0, 0.0],"T":[0.0, 0.0, 1.0, 0.0],"G":[0.0, 0.0, 0.0, 1.0]}
    for aBase in seqIn:
        encSeq+=seqToOneHotDict[aBase]  
    encSeq=np.array(encSeq)
    return encSeq   

# def onehot_conversion_sequence(cls, letter):
#         one_hot_map = {
#             "A": np.asarray([1, 0, 0, 0],dtype=np.float32), "a": np.asarray([1, 0, 0, 0],dtype=np.float32),
#             "C": np.asarray([0, 1, 0, 0],dtype=np.float32), "c": np.asarray([0, 1, 0, 0],dtype=np.float32),
#             "G": np.asarray([0, 0, 1, 0],dtype=np.float32), "g": np.asarray([0, 0, 1, 0],dtype=np.float32),
#             "T": np.asarray([0, 0, 0, 1],dtype=np.float32), "t": np.asarray([0, 0, 0, 1],dtype=np.float32),
#             "N": np.asarray([0, 0, 0, 0],dtype=np.float32), "n": np.asarray([0, 0, 0, 0],dtype=np.float32)}
#         return one_hot_map[letter]

# def one_hot_encoder(self):
#     tmp = []
#     for letter in self.sequence:
#         tmp.append(self.onehot_conversion_sequence(letter))
#     out = np.vstack(tmp)
#     return (out)

def dnaCheck(sequence):
    dna = set('ACTG')
    return all(base.upper() in dna for base in sequence)