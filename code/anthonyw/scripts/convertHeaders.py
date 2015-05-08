#!/usr/bin/python3
#
# Convert FASTQ paired headers to include both alignments


import os
import subprocess
import argparse
import re
from concurrent.futures.thread import ThreadPoolExecutor

def fixHeader(passPrefix):
    # Open FASTQs for reading
    fastqName1 = "{0}-1.fq".format(passPrefix)
    fastqName2 = "{0}-2.fq".format(passPrefix)
    fastqHandle1 = open(fastqName1, 'r')
    fastqHandle2 = open(fastqName2, 'r')
    
    # Open temporary file for writing (with .tmp extension)
    newName1 = "1.tmp"
    newName2 = "2.tmp"
    newHandle1 = open(newName1, 'w')
    newHandle2 = open(newName2, 'w')
        
    # Iterate through each line in FASTQ (without reading into memory)
    sectionIdx = 0
    for fastqLine1 in fastqHandle1:
        # Readcorresponding line for file 2
        fastqLine2 = fastqHandle2.readline()
        
        # For header, generate and write new sequence header
        if sectionIdx == 0:
            headerMatch1 = re.match('^(.*):-:([0-9]+)/1$', fastqLine1)
            headerMatch2 = re.match('^(.*):-:([0-9]+)/2$', fastqLine2)
            
            newHandle1.write("{0}:-:{1}:-:{2}/1\n".format(headerMatch1.group(1), headerMatch1.group(2), headerMatch2.group(2)))
            newHandle2.write("{0}:-:{1}:-:{2}/2\n".format(headerMatch2.group(1), headerMatch1.group(2), headerMatch2.group(2)))
                        
        else:
            # Not a header, simply copy line
            newHandle1.write(fastqLine1)
            newHandle2.write(fastqLine2)
            
        # Advance additional SAM entry for PE files, and FASTQ section index
        sectionIdx = (sectionIdx + 1) % 4

    # Close files
    fastqHandle1.close()
    fastqHandle2.close()
    newHandle1.close()
    newHandle2.close()
    
    # Delete original FASTQ, replace with temporary file
    os.remove(fastqName1)
    os.remove(fastqName2)
    os.rename(newName1, fastqName1)
    os.rename(newName2, fastqName2)


fastqRoot, fastqDirs, fastqFiles = next(os.walk(os.curdir))
        
for currentFile in fastqFiles:
    if "-1.fq" in currentFile:
       print("Converting {0}...".format(currentFile))
       fixHeader(os.path.splitext(currentFile)[0][:-2]) 
