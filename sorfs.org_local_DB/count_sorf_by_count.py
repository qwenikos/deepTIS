initDir="../"
import os

dbFileName=initDir+"data/sorfs.org/sorfs_allData2022.tab"


workingDir="orfDBWithTISFlanking/"
command="mkdir -p "+workingDir
os.system(command)



    
#####################################

dbFile=open(dbFileName,"r")
header=dbFile.readline()
headerCols= header.rstrip().split("\t")
colsIdx=0
for aHeader in headerCols:
    
    print (str(colsIdx)+") "+str(aHeader))
    colsIdx+=1

aTempList=[]
allCnt=0
mitochondrianCnt=0
selectedSizeCnt=0
minLength=18
maxLength=100
aSorfIdDict={}


aSorfAnnotationList=["sORF"]
sorfAnnotatedCnt=0
intronicAnnotatedCnt=0
uknownSizeCnt=0
withNoMetricsCnt=0
sorfAnnotatedWithMetricsCnt=0
intronicAnnotatedWirthMetricsCnt=0
for aSorf in dbFile:
    allCnt+=1
    cols=aSorf.rstrip().split("\t")
    aSorfId=cols[0]
    if aSorfId not in aSorfIdDict:
        aSorfIdDict[aSorfId]=""
    aSorfChrom=cols[2]
    aSorfChrom="chr"+aSorfChrom
    if aSorfChrom=="chrMT":
        mitochondrianCnt+=1
        continue
    aSorfStrand=cols[3]
    if (aSorfStrand=="1"):
        aSorfStrand="+"
    elif (aSorfStrand=="-1"):
        aSorfStrand="-"
        
    aSorfStart=int(cols[4])
    aSorfEnd=int(cols[5])
    aSorfLength=int(cols[9])
    aSorfTransSeq=cols[13]
    aSorfFlossClass=cols[28]
    aSorfStartCodon=cols[10]
    aSorfPcExonsOverlap=cols[24]
    aSorfEnsemblId=cols[25]
    aSorfFLOSSScore=cols[27]
    aSorfInFrame=cols[23]
    aSorfAnnotation=cols[16]
    aSorfFlosslId=cols[27]
    aSorfOrfscore=cols[29]
    aSorfPhastCon=cols[31]
    aSorfPhyloP=cols[32]
    sorfStart=aSorfStart-1 ##zero based to 1 based
    sorfEnd=aSorfEnd

    if ((aSorfLength>=minLength) and (aSorfLength<=maxLength)):
        selectedSizeCnt+=1
    else:
        uknownSizeCnt+=1
        continue

    if aSorfAnnotation=="sORF":
        sorfAnnotatedCnt+=1

    if aSorfAnnotation=="intronic":
        intronicAnnotatedCnt+=1

    if aSorfOrfscore=="" or aSorfPhastCon=="" or  aSorfPhyloP=="" or aSorfFlosslId=="":
        withNoMetricsCnt+=1
    else:
        if aSorfAnnotation=="sORF":
            sorfAnnotatedWithMetricsCnt+=1
        if aSorfAnnotation=="intronic":
            intronicAnnotatedWirthMetricsCnt+=1


print ("__________________________")
print ("allCnt",allCnt)
print ("selectedSizeCnt",selectedSizeCnt,minLength,maxLength)
print ("mitochondrianCnt",mitochondrianCnt)
print ("len(aSorfIdDict)",len(aSorfIdDict))
print ("sorfAnnotatedCnt",sorfAnnotatedCnt)
print ("uknownSizeCnt",uknownSizeCnt)
print ("intronicAnnotatedCnt",intronicAnnotatedCnt)
print ("withNoMetricsCnt",withNoMetricsCnt)
print ("sorfAnnotatedWithMetricsCnt",sorfAnnotatedWithMetricsCnt)
print ("intronicAnnotatedWirthMetricsCnt",intronicAnnotatedWirthMetricsCnt)