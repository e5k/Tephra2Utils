#!/usr/bin/python
from string import Template
from sys import argv
import re

# Read parameters.README
def splt(line):
    ln = line.split(':')
    ln = ln[1].split()
    return ln[0]

file = open("parameters.README", "r")
for line in file:
    if re.search("Max Column Height", line):
        ht = splt(line)
    elif re.search("Total Mass Ejected:", line):
        mass = splt(line)
    elif re.search("Alpha Param:", line):
        alpha = splt(line)
    elif re.search("Beta Param:", line):
        beta = splt(line)
    elif re.search("Diffusion Coefficient:", line):
        diff = splt(line)
    elif re.search("Fall Time Threshold:", line):
        FTT = splt(line)
    elif re.search("Eddy Constant:", line):
        eddy = splt(line)
    elif re.search("Max Particle Size:", line):
        minPhi = splt(line)
    elif re.search("Min Particle Size:", line):
        maxPhi = splt(line)
    elif re.search("Median Size:", line):
        medPhi = splt(line)
    elif re.search("Std. Dev.:", line):
        sigPhi = splt(line)
    elif re.match("Parameter Ranges:", line):
        break

# Read inversionConfig.conf
def splt2(line):
    ln = line.split('=')
    ln = ln[1].split('\n')
    return ln[0]

file = open("../inversionConfig.conf", "r")
for line in file:
    if re.search("ventE", line):
        ventE = splt2(line)
    elif re.search("ventN", line):
        ventN = splt2(line)
    elif re.search("ventA", line):
        ventA = splt2(line)
    elif re.search("lithicDensity", line):
        lithicDensity = splt2(line)
    elif re.search("pumiceDensity", line):
        pumiceDensity = splt2(line)
    elif re.search("colSteps", line):
        colSteps = splt2(line)
    elif re.search("partSteps", line):
        partSteps = splt2(line)
    elif re.search("plumeModel", line):
        plumeModel = splt2(line)

# Write file
with open( '../../_templates/forwardConfTemplate.conf', 'r' ) as inF:
  config = Template(inF.read() )
  conf = config.substitute(ventE=ventE, ventN=ventN, ventA=ventA, ht=ht, mass=mass, alpha=alpha, beta=beta, minPhi=minPhi, maxPhi=maxPhi, medPhi=medPhi, sigPhi=sigPhi, eddy=eddy, diff=diff, FTT=FTT, lithicDensity=lithicDensity, pumiceDensity=pumiceDensity, colSteps=colSteps, partSteps=partSteps, plumeModel=plumeModel)
  with open( 'tephra2.conf', 'w' ) as outF:
    outF.write( conf )



