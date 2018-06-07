#!/usr/bin/python

from string import Template
from sys import argv

minHt = argv[1]
maxHt = argv[2]
minMass = argv[3]
maxMass = argv[4] 
ventE = argv[5] 
ventN = argv[6] 
ventA = argv[7] 
minDiff = argv[8] 
maxDiff = argv[9] 
eddy = argv[10]
minMedPhi = argv[11]
maxMedPhi = argv[12]
minSigPhi = argv[13]
maxSigPhi = argv[14]
minAlpha = argv[15]
maxAlpha = argv[16]
minBeta = argv[17]
maxBeta = argv[18]
minFTT = argv[19]
maxFTT = argv[20]
minWindSpeed = argv[21]
maxWindSpeed = argv[22]
minWindDir = argv[23]
maxWindDir = argv[24]
plumeModel = argv[25]
fixedWind = argv[26]
windLevels = argv[27]
colSteps = argv[28]
partSteps = argv[29]
lithicDensity = argv[30]
pumiceDensity = argv[31]
minPhi = argv[32]
maxPhi = argv[33]
fitTest = argv[34]

with open( '../../_templates/inversionConfTemplate.conf', 'r' ) as inF:
  config = Template(inF.read() )
  conf = config.substitute(minHt=minHt, maxHt=maxHt, minMass=minMass, maxMass=maxMass, ventE=ventE, ventN=ventN, ventA=ventA, minDiff=minDiff, maxDiff=maxDiff, eddy=eddy, minMedPhi=minMedPhi, maxMedPhi=maxMedPhi, minSigPhi=minSigPhi, maxSigPhi=maxSigPhi, minAlpha=minAlpha, maxAlpha=maxAlpha, minBeta=minBeta, maxBeta=maxBeta, minFTT=minFTT, maxFTT=maxFTT, minWindSpeed=minWindSpeed, maxWindSpeed=maxWindSpeed, minWindDir=minWindDir, maxWindDir=maxWindDir, plumeModel=plumeModel, fixedWind=fixedWind, windLevels=windLevels, colSteps=colSteps, partSteps=partSteps, lithicDensity=lithicDensity, pumiceDensity=pumiceDensity, minPhi=minPhi, maxPhi=maxPhi, fitTest=fitTest)
  with open( 'tmp.conf', 'w' ) as outF:
    outF.write( conf )