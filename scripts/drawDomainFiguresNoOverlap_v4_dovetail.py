import sys, re, os
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import random

sys.path.append('/nfs0/BB/Hendrix_Lab/Hops/svgwrite/svgwrite-1.2.1/')
import svgwrite
from svgwrite import cm, mm

def readGFF(geneGFF):
    strandInfo = {}
    with open(geneGFF,'r') as GFF:
        for line in GFF:
            if not line.startswith('#'):
                contigID,source,feature,start,stop,score,strand,frame,attribute  = line.strip().split("\t")
                if feature == 'mRNA':
                    getGeneID = re.search('ID=(.+);Parent',attribute)
                    geneID = getGeneID.group(1)
		    #print geneID
                    strandInfo[geneID] = strand
    return(strandInfo)


def readHMMFile(hmmFile,eValueThresh,scaleFactor):
    hmmCoordDict = {}
    with open(hmmFile,'r') as HMM:
        for line in HMM:
            if not line.startswith('#'):
                if not line.startswith('-'):
                    if not line.isspace():
                        domainInfo = line.strip().split()[0:22]
                        accessionID = domainInfo[1]
                        hopGeneID = domainInfo[3]
                        hopGeneLen = domainInfo[5]
                        eValue = domainInfo[6]
                        score = domainInfo[7]
                        domainStart = domainInfo[17]
                        domainStop = domainInfo[18]
                        alignmentProbability = domainInfo[21]
                        accessionDesc = line.strip().split()[22:]
                        accessionDesc = ' '.join(accessionDesc)
                        eValue = float(eValue)
                        hopGeneLen = int(hopGeneLen)
                        scaledHopGeneLen = float(hopGeneLen) / scaleFactor
                        domainStart = int(domainStart)
                        domainStop = int(domainStop)
                        scaledStart = float(domainStart) / scaleFactor
                        scaledStop = float(domainStop) / scaleFactor
                        if eValue < eValueThresh:
                            if (hopGeneID,hopGeneLen,scaledHopGeneLen) not in hmmCoordDict:
                                hmmCoordDict[(hopGeneID,hopGeneLen,scaledHopGeneLen)] = []
                            hmmCoordDict[(hopGeneID,hopGeneLen,scaledHopGeneLen)].append((accessionID,accessionDesc,int(domainStart),int(domainStop),int(scaledStart),int(scaledStop)))
    return(hmmCoordDict)

def getPosCount(countArray,start,stop,updatedMaxCount):
    for i in range(start,stop):
        countArray[i] = updatedMaxCount
    return(countArray)

def getMaxCountValue(countArray,start,stop):
    maxValue = 0
    for i in range(start,stop):
        if maxValue < countArray[i]:
            maxValue = countArray[i]
    return(maxValue)

def assessOverlap(hmmCoordDict):
    filteredGeneDict = {}
    for hopGeneID,hopGeneLen,scaledHopGeneLen in hmmCoordDict:
        maxCountList = []
        scaledGeneLen = int(scaledHopGeneLen)
        # countArray = [0]*scaledGeneLen
        countArray = [0]*hopGeneLen
        
        for accessionID,accessionDesc,domainStart,domainStop,scaledStart,scaledStop in hmmCoordDict[(hopGeneID,hopGeneLen,scaledHopGeneLen)]:
            maxCount = overlap(countArray,domainStart,domainStop,maxCountList)
            maxCountList.append(maxCount)
        maximumCount = max(maxCountList)
        if maximumCount == 1:
            # print(hopGeneID,maximumCount)
            filteredGeneDict[hopGeneID] = maximumCount
    return(filteredGeneDict)

def countDomains(hmmCoordDict,filteredGeneDict):
    domainCountDict = {}
    for hopGeneID,hopGeneLen,scaledHopGeneLen in hmmCoordDict:
        if hopGeneID in filteredGeneDict:
            accessionIDDict = {}
            domainCount = 0
            for accessionID,accessionDesc,domainStart,domainStop,scaledStart,scaledStop in hmmCoordDict[(hopGeneID,hopGeneLen,scaledHopGeneLen)]:
                #print(hopGeneID,accessionID)
                if accessionID not in accessionIDDict:
                    accessionIDDict[accessionID] = 1
                    domainCount += 1
                    domainCountDict[hopGeneID] = domainCount
    return(domainCountDict)

def overlap(countArray,domStart,domStop,maxCountList):
    maxCount = getMaxCountValue(countArray,domStart,domStop)
    updatedMaxCount = maxCount + 1
    countArray = getPosCount(countArray,domStart,domStop,updatedMaxCount)
    maxCount = getMaxCountValue(countArray,domStart,domStop)
    return(maxCount)

def createColorPalette(geneDict,geneID,domainCount):
    colors = ['#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061']
    colorStuff = {}
    accessionDict = {}
    accessionList = []
    for accessionID,accessionDesc,domainStart,domainStop,scaledStart,scaledStop in geneDict:
        if accessionID not in accessionDict:
            accessionDict[accessionID] = 1
            accessionList.append(accessionID)
    for i in range(len(accessionList)):
        if accessionList[i] not in colorStuff:
            colorStuff[accessionList[i]] = colors[i]
            #print(accessionList[i],colors[i])
    return(colorStuff)

def drawGene(hmmCoordDict,filteredGeneDict,domainCountDict,strandInfo):
    for hopGeneID,hopGeneLen,scaledHopGeneLen in hmmCoordDict:
        #hmmCoordDict[(hopGeneID,hopGeneLen,scaledHopGeneLen)].sort(key=lambda x:x[2], reverse=False)
        if hopGeneID in filteredGeneDict and hopGeneID in domainCountDict:
            domCount = domainCountDict[hopGeneID]
            colorDict = {}
            accessionDict = {}
            otherAccessionDict = {}
            domainState = 0
            legendTextOffset = 0
            colorCount = 0
            scaledGeneLen = int(scaledHopGeneLen)
            figWidth = 1800
            figHeight = 35*domCount + 100
            fig = svgwrite.Drawing(filename=hopGeneID + ".svg", size=(figWidth, figHeight), profile='full')
            # countArray = [0]*scaledGeneLen
            countArray = [0]*hopGeneLen

            geneDrawingLength = 1000
            #fig.add(fig.rect((1,52.5), (scaledGeneLen,30),fill='#2166ac', rx=2, ry=2))
            fig.add(fig.rect((1,52.5), (geneDrawingLength,30),fill='#2166ac', rx=2, ry=2))
            # fig.add(fig.rect((0, (100)), (30,12),fill='#2166ac', rx=1, ry=1))
            fig.add(fig.rect((0, (100)), (30,12),fill='#2166ac', rx=1, ry=1))
            fig.add(fig.text(hopGeneID,
                             insert=(40, (111.75)),
                             fill='black', font_size='16px')
                    )
            
            accessions = {}
            printDict = {}
            OUT = open(hopGeneID + '.pfam.txt','w')
            if domCount <= 10:
                colorStuff = createColorPalette(hmmCoordDict[(hopGeneID,hopGeneLen,scaledHopGeneLen)],hopGeneID,domCount)
                for accessionID,accessionDesc,domainStart,domainStop,scaledStart,scaledStop in hmmCoordDict[(hopGeneID,hopGeneLen,scaledHopGeneLen)]:
                    # print(hopGeneID,accessionID,accessionDesc,domainStart,domainStop,scaledStart,scaledStop)
                    if accessionID not in printDict:
                        printDict[accessionID] = 1
                        #print hopGeneID,accessionID,accessionDesc
                        OUT.write("%s\t%s\t%s\n" % (hopGeneID,accessionID,accessionDesc))
                    if accessionID not in accessionDict:
                        domainState += 1.15
                        accessionDict[accessionID] = domainState
    
                    if accessionID in colorStuff:
                        domainColor = colorStuff[accessionID]
                        accessionInfo = accessionID + "\t" + accessionDesc

                        if accessionID not in accessions:
                            accessions[accessionID] = 1
                            legendTextOffset += 17
                            colorState = legendTextOffset
                            
                            legendOffset = scaledGeneLen + 25
                            
                            fig.add(fig.rect((0, (100+colorState)), (30,12),fill=domainColor, rx=1, ry=1))
                            fig.add(fig.text(accessionInfo,
                                             insert=(40, 111.75+colorState),
                                             fill='black', font_size='16px')
                                    )
                        
                    drawDomains(hopGeneID,countArray,domainStart,domainStop,hopGeneLen,accessionID,accessionDesc,domainColor,fig,accessionInfo,otherAccessionDict,figHeight,geneDrawingLength,strandInfo)
            if domCount > 10:
                for accessionID,accessionDesc,domainStart,domainStop,scaledStart,scaledStop in hmmCoordDict[(hopGeneID,hopGeneLen,scaledHopGeneLen)]:
                    # print(hopGeneID,accessionID,accessionDesc,domainStart,domainStop,scaledStart,scaledStop)
                    if accessionID not in printDict:
                        printDict[accessionID] = 1
                        #print hopGeneID,accessionID,accessionDesc
                        OUT.write("%s\t%s\t%s\n" % (hopGeneID,accessionID,accessionDesc))
                    if accessionID not in accessionDict:
                        domainState += 1.15
                        accessionDict[accessionID] = domainState
                        
                    if accessionID not in colorDict:
                        accessionInfo = accessionID + "\t" + accessionDesc
                        r = lambda: random.randint(0,255)
                        randomColor = ('#%02X%02X%02X' % (r(),r(),r()))
                        accessionInfo = accessionID + "\t" + accessionDesc
                        colorDict[accessionID] = randomColor
                        legendTextOffset += 17
                        colorState = legendTextOffset
                        legendOffset = scaledGeneLen + 25
                        
                        fig.add(fig.rect((0, (100+colorState)), (30,12),fill=randomColor, rx=1, ry=1))
                        fig.add(fig.text(accessionInfo,
                                         insert=(40, 111.75+colorState),
                                         fill='black', font_size='16px')
                                )
                    drawDomains(hopGeneID,countArray,domainStart,domainStop,hopGeneLen,accessionID,accessionDesc,randomColor,fig,accessionInfo,otherAccessionDict,figHeight,geneDrawingLength,strandInfo)


def drawDomains(geneID,countArray,domStart,domStop,geneLen,accID,accDesc,domainColor,fig,accessionInfo,otherAccessionDict,figHeight,geneDrawingLength,strandInfo):
    alignLen = abs(domStop-domStart + 1)
    maxCount = getMaxCountValue(countArray,domStart,domStop)
    shapes = fig.add(fig.g(id='shapes', fill='darkmagenta'))
    strand = strandInfo[geneID]

    domainDrawingStart = float(geneDrawingLength*domStart)/geneLen
    domainDrawingLength = float(geneDrawingLength*alignLen)/geneLen
    domainDrawingStopPos = domainDrawingStart + domainDrawingLength
    #print domainDrawingStart,domainDrawingStopPos,domainDrawingLength

    # fig.add(fig.rect((domStart,(maxCount*(-35))-2.5), (alignLen,35),fill=domainColor,rx=2, ry=2))
    ## fig.add(fig.rect((domStart,(maxCount*(-35))+50), (alignLen,35),fill=domainColor,rx=2, ry=2))
    if strand == '-':
        revStrandStart = geneDrawingLength - domainDrawingStopPos
        fig.add(fig.rect((revStrandStart,(maxCount*(-35))+50), (domainDrawingLength,35),fill=domainColor,rx=2, ry=2))

        updatedMaxCount = maxCount + 1
        countArray = getPosCount(countArray,domStart,domStop,updatedMaxCount)
        maxCount = getMaxCountValue(countArray,domStart,domStop)
    
    else:
        fig.add(fig.rect((domainDrawingStart,(maxCount*(-35))+50), (domainDrawingLength,35),fill=domainColor,rx=2, ry=2))

        updatedMaxCount = maxCount + 1
        countArray = getPosCount(countArray,domStart,domStop,updatedMaxCount)
        maxCount = getMaxCountValue(countArray,domStart,domStop)

    fig.save()


########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <hmmscan domtblout file> <gene gff> \n"
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

hmmFile = sys.argv[1]
geneGFF = sys.argv[2]

eValueThresh = 1e-5
scaleFactor = 10

strandInfo = readGFF(geneGFF)
hmmCoordDict = readHMMFile(hmmFile,eValueThresh,scaleFactor)
filteredGeneDict = assessOverlap(hmmCoordDict)
domainCountDict = countDomains(hmmCoordDict,filteredGeneDict)

drawGene(hmmCoordDict,filteredGeneDict,domainCountDict,strandInfo)
