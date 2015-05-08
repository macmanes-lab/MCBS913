#!/usr/bin/python3
#
# Convert all SAM statistics in a directory to a CSV file

import os
import csv
import re
    
# Open CSV file for writing
with open('samstats.csv', 'w') as csvHandle:
    csvWriter = csv.writer(csvHandle, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    
    # Write header
    csvWriter.writerow(['Filename', 'Total Reads (QC-P)', 'Total Reads (QC-F)', 'Mapped (QC-P)', 'Mapped (QC-F)', 'Mapped % (QC-P)'])
    
    # Write data from each file
    samRoot, samDirs, samFiles = next(os.walk(os.curdir))
            
    for currentFile in samFiles:
        if '.samstat' in currentFile:
            samHandle = open(currentFile, 'r')
    
            rowData = list()
            
            # Filename
            rowData.append(os.path.basename(currentFile))
            
            # Total Reads (Success and Failure)
            samLine = samHandle.readline()
            samMatch = re.match('^([0-9]+) \\+ ([0-9]+).*$', samLine)
            rowData.append(samMatch.group(1))
            rowData.append(samMatch.group(2))
            
            # Mapped (Success, Failure, Percentage)
            samLine = samHandle.readline()
            samLine = samHandle.readline()
            samMatch = re.match('^([0-9]+) \\+ ([0-9]+).*\\(([0-9]+\\.[0-9]+).*\\)$', samLine)
            rowData.append(samMatch.group(1))
            rowData.append(samMatch.group(2))
            rowData.append(samMatch.group(3))
            
            # Write Row
            csvWriter.writerow(rowData)
