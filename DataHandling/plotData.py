import os
import sys
import matplotlib.pyplot as plot
from parseAsc import *

xLabel = "Number of Molecules"
yLabel = "Time (seconds)"
title = "Non-Bonded Molecule Simulation"

def readFile(fileName):
    '''Returns a string containing the contents of a file.'''
    contents = open(fileName).read()
    print contents
    return contents

def getSize(fileName):
    '''Gets the size of the sim from the filename.'''
    size = '0'
    fileName = fileName.split('/')
    fileName = fileName[1]
    print fileName
    for i in fileName:
        if i in string.digits:
            size += i
            print size
        else:
            break
    return int(size)

def plotData(linear, parallel):
    parallel[500] = 68789
    print parallel
    print linear
    pX = []
    pY = []

    lX = []
    lY = []

    for i in linear :
        lX.append(i)
        lY.append(linear[i])
    for i in parallel:
        pX.append(i)
        pY.append(parallel[i])
    lX.sort()
    lY.sort()
    pX.sort()
    pY.sort()

    print lX
    print lY
    print pX
    print pY

    plot.title(title)
    plot.xlabel(xLabel)
    plot.ylabel(yLabel)
    plot.loglog(lX, lY, color='green', linestyle='solid', marker='o', label='Linear')
    plot.loglog(pX, pY, color='blue', linestyle='solid', marker='o', label='Parallel')
    plot.legend(loc='upper left')
    plot.savefig("cycle2Data.png")

def main():
    linearDir = sys.argv[1]
    parallelDir = sys.argv[2]
    
    linearFiles = os.listdir(linearDir)
    parallelFiles = os.listdir(parallelDir)
    
    linearTime = {}
    parallelTime = {}

    for i in linearFiles:
        data = readFile(linearDir + i)
        size = getSize(linearDir + i)
        cpuTime = getCPUTime(data)
        linearTime[size] = cpuTime

    for i in parallelFiles:
        data = readFile(parallelDir + i)
        size = getSize(parallelDir + i)
        cpuTime = getCPUTime(data)
        parallelTime[size] = cpuTime

    print linearTime
    print parallelTime

    plotData(linearTime, parallelTime)

if __name__ == '__main__':
    main()
