#!/usr/bin/env python
__author__ = 'Guanwen Wang'

# Guanwen Not part of the final project
import sys, os, fnmatch

# get Assembly for the single and pair ended genomes

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    MAGENTA = "\033[.35m"
    ENDCODE = "\033[.0m"

PrintAll = 0

TotalFound = 0
result = []

HHListS = 'HHListS.txt'
HHListP = 'HHListSingleP.txt'

pathArr = ['/thomas1/mcbs913/shared_data1/datasets/reads']
patternS = 'Ha*-s*10.*'
#patternP = 'Ha*-*-*-*-10-*-*'
patternP = 'Ha*-*-*-*-10-*-1.*'

def GenomeFilter(pattern, pathArr):
    global result 
    global TotalFound
    for path in pathArr:
        for root, dirs, files in os.walk(path):
            found = False
            for name in files:
                if fnmatch.fnmatch(name, pattern):
                    found = True
                    targetFile = os.path.join(root, name)
                    if PrintAll == 1:
                        print targetFile
		    TotalFound += 1
                    result.append(targetFile)
            if( found == False):
                print root 
                print bcolors.FAIL + "WARNING: No Genome files found\n" + bcolors.ENDC
   


GenomeFilter(patternS, pathArr)
result.sort()

with open(HHListS, 'w') as f:
    for s in result:
        f.write(s + '\n');

if PrintAll == 1:
   print TotalFound
   print result

TotalFound = 0
del result[:]

GenomeFilter(patternP, pathArr)
result.sort()

with open(HHListP, 'w') as f:
    for s in result:
        f.write(s + '\n');

if PrintAll == 1:
   print TotalFound
   print result
