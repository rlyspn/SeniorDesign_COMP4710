import re
import string
def getCPUTime(fileStr):
    '''Gets the cpu time (sec) from an asc output file.'''
    timeReg = re.compile('cput=(\d\d):(\d\d):(\d\d),')
    usedTime = re.findall(timeReg, fileStr)[1]
    return int(usedTime[0]) * 3600 + int(usedTime[1]) * 60 + int(usedTime[2])

def getWallTime(fileStr):
    '''Gets the wall time (sec) from an asc output file.'''
    timeReg = re.compile('walltime=(\d\d):(\d\d):(\d\d)$')
    usedTime = re.findall(timeReg, fileStr)[1]
    return int(usedTime[0]) * 3600 + int(usedTime[1]) * 60 + int(usedTime[2])

def getMemory(fileStr):
    '''Gets the memory (kb) usage from an asc output file.'''
    memReg = re.compile(',mem=(\d+)kb,')
    mem = re.findall(memReg, fileStr)
    return int(mem[0])
    
def getVirtualMemory(fileStr):
    '''Returns the amount of virtual memory used in kilobytes.'''
    vMemReg = re.complie(',vmem=(\d+)kb,')
    vMem = re.findall(vMemReg, fileStr)
    return int(vMem[0])


