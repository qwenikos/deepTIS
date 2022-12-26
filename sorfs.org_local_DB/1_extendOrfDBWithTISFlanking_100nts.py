initDir="../"
import os

dbFileName=initDir+"data/sorfs.org/sorfs_allData2022.tab"


workingDir="orfDBWithTISFlanking_100nts/"
command="mkdir -p "+workingDir
os.system(command)

sOrfsAllFieldsWithDublicatesBedFileName        =workingDir+"1_sOrfsAllFields_WithDublicates.bed"
sOrfsSelectedFieldsWithDublicatesBedFileName   =workingDir+"2_sOrfsSelectedFields_WithDublicates.bed"

sOrfsAllFieldsBedFileName            =workingDir+"3_sOrfsAllFields.bed"
sOrfsSelectedFieldsBedFileName       =workingDir+"4_sOrfsSelectedFields.bed"

sOrfsLengthFileName                     =workingDir+"5_sOrfsLength.csv"
sOrfsWithExtendedFlankBedFileName       =workingDir+"sOrfsWithExtendedFlank.bed" ##with TIS flank
sOrfsWithExtendedFlankTabFileName       =workingDir+"sOrfsWithExtendedFlank.tab"
sOrfsWithExtendedFlankFaFileName        =workingDir+"sOrfsWithExtendedFlank.fa"

sOrfs_sORF_BedFileName      =workingDir+"sOrfs_sORF.bed"
sOrfs_intronic_BedFileName  =workingDir+"sOrfs_intronic.bed"
sOrfs_sORF_FaFileName      =workingDir+"sOrfs_sORF.fa"
sOrfs_intronic_FaFileName  =workingDir+"sOrfs_intronic.fa"
# outputDBFileName                        =workingDir+"sOrfsDBWithTIS_100nts.tab"
flankingInterval=100
seqDict={}
dublicatesCnt=0
#####################################

refGenome="hg38"
genomeDir=initDir+"data/GENOME/"

if refGenome=="hg38":
    genomeFileName=genomeDir+"hg38.fa"
    chromSizes=genomeDir+"hg38.chrom.sizes"
    
if refGenome=="hg19":
    genomeFileName="/mnt/raid0/raid_common/Genomes/hg19/hg19.fa"
    chromSizes="/mnt/raid0/raid_common/Genomes/hg19/chrom.sizes"
    
#####################################
######create the bed file

dbFile=open(dbFileName,"r")
header=dbFile.readline()
headerCols= header.rstrip().split("\t")
colsIdx=0
for aHeader in headerCols:
    
    print (str(colsIdx)+") "+str(aHeader))
    colsIdx+=1

sOrfsAllFieldsWithDublicatesBedFile=open(sOrfsAllFieldsWithDublicatesBedFileName,"w")
sOrfsSelectedFieldsWithDublicatesBedFile=open(sOrfsSelectedFieldsWithDublicatesBedFileName,"w")
sOrfsAllFieldsBedFile=open(sOrfsAllFieldsBedFileName,"w")
sOrfsSelectedFieldsBedFile=open(sOrfsSelectedFieldsBedFileName,"w")
sOrfsLengthFile=open(sOrfsLengthFileName,"w")

sOrfs_sORF_BedFile=open(sOrfs_sORF_BedFileName,"w")
sOrfs_intronic_BedFile=open(sOrfs_intronic_BedFileName,"w")

aTempList=[]
cnt1=0
cnt2=0
cnt3=0
for aSorf in dbFile:
    cnt1+=1
    cols=aSorf.rstrip().split("\t")
    aSorfId=cols[0]
    aSorfChrom=cols[2]
    aSorfChrom="chr"+aSorfChrom
    if aSorfChrom=="chrMT":
        cnt3+=1
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
    #print (aSorfInFrame)
    aSorfStartCodonList=["ATG"]
    aSorfFlossClassList=["Good"]
    aSorfFlossClassList=["Not in cutoff range"]
    aSorfFlossClassList=["No reads"]
    aSorfFlossClassList=["Good","Extreme"]
    aSorfPcExonsOverlapList=["Yes","No",""]
    aSorfEnsemblIdList=[""]
    aSorfInFrameList =["Yes","No",""]
    aSorfAnnotationList=["sORF"]
    sorfStart=aSorfStart-1 ##zero based to 1 based
    sorfEnd=aSorfEnd
    
    newName="|".join(map(str,cols))
    selectedFieldsName="|".join([aSorfId,aSorfAnnotation])
    #print (newName)
    aAllBedLine="\t".join(map(str,[aSorfChrom,aSorfStart-1,aSorfEnd,newName,aSorfFLOSSScore,aSorfStrand]))
    aSelectedBedLine="\t".join(map(str,[aSorfChrom,aSorfStart-1,aSorfEnd,selectedFieldsName,"0",aSorfStrand]))
    if aSorfEnd>aSorfStart:
        sOrfsAllFieldsWithDublicatesBedFile.write(aAllBedLine+"\n")
        sOrfsSelectedFieldsWithDublicatesBedFile.write(aSelectedBedLine+"\n")
        cnt2+=1
    else:
        print (">>",aSorfStart,"^",aSorfEnd,"<<")

    
    if not (aSorfTransSeq in seqDict):
        seqDict[aSorfTransSeq]=[aSorfId]
        sOrfsAllFieldsBedFile.write(aAllBedLine+"\n")
        sOrfsSelectedFieldsBedFile.write(aSelectedBedLine+"\n")
        sOrfsLengthFile.write(str(aSorfLength)+"\n")
        if (aSorfAnnotation=="sORF"):
            sOrfs_sORF_BedFile.write(aAllBedLine+"\n")
        if (aSorfAnnotation=="intronic"):
            sOrfs_intronic_BedFile.write(aAllBedLine+"\n")

    else:
        dublicatesCnt+=1
        seqDict[aSorfTransSeq]+=[aSorfId]

print (cnt1,cnt2,cnt3)
sOrfsAllFieldsWithDublicatesBedFile.close()
sOrfsSelectedFieldsWithDublicatesBedFile.close()
sOrfsAllFieldsBedFile.close()
sOrfsSelectedFieldsBedFile.close()
sOrfsLengthFile.close()
sOrfs_sORF_BedFile.close()
sOrfs_intronic_BedFile.close()
#############################################

command="bedtools slop  -s -g "+chromSizes+" -i "+sOrfsSelectedFieldsBedFileName+"  -l "+str(flankingInterval)+" -r 0 > "+sOrfsWithExtendedFlankBedFileName;
print (command)
os.system(command)

#############################################

command="bedtools getfasta -name -tab -s -fi "+genomeFileName+" -bed "+sOrfsWithExtendedFlankBedFileName+" -fo "+sOrfsWithExtendedFlankTabFileName;
print (command)
os.system(command)

#############################################

command="bedtools getfasta -name  -s -fi "+genomeFileName+" -bed "+sOrfsWithExtendedFlankBedFileName+" -fo "+sOrfsWithExtendedFlankFaFileName;
print (command)
os.system(command)

#############################################

command="bedtools getfasta -name  -s -fi "+genomeFileName+" -bed "+sOrfs_sORF_BedFileName+" -fo "+sOrfs_sORF_FaFileName;
print (command)
os.system(command)

#############################################

command="bedtools getfasta -name  -s -fi "+genomeFileName+" -bed "+sOrfs_intronic_BedFileName+" -fo "+sOrfs_intronic_FaFileName;
print (command)
os.system(command)

#############################################

exit()
##reconstruct DB with extra columns 34)TIS and 35)flanking interval length
# sOrfWithExtendedFlankTabFile=open(sOrfsWithExtendedFlankTabFileName,"r")
# outputDBFile=open(outputDBFileName,"w")
# newHedaderCols=headerCols+["TisSequence","StartPos"]
# print (newHedaderCols)

# outputDBFile.write("\t".join(newHedaderCols)+"\n")
# for aLine in sOrfWithExtendedFlankTabFile:

#     cols=aLine.rstrip().split("\t")
#     names=cols[0]
#     tisSequence=cols[1]
#     tisSequence=tisSequence[flankingInterval-7:flankingInterval+5]
#     originalCols=names.split("|")
#     newColsList=originalCols+[tisSequence.upper()]+[str(8)]
#     newLine="\t".join(map(str,newColsList))
#     outputDBFile.write(newLine+"\n")
#     #print (newLine)
    
# sOrfWithExtendedFlankTabFile.close()
# outputDBFile.close()



    
