#!/usr/bin/python3
#
# Read Simulator via ART - readSimulator.py
#
# Anthony Westbrook             
# Seth Hager  
#
# usage: readSimulator.py [-h] [--art ART] [--samtools SAMTOOLS] [--gzip GZIP]
#                         [--fastadir FASTADIR] [--fasta FASTA]
#                         [--output OUTPUT] [--lengths LENGTHS]
#                         [--fragments FRAGMENTS] [--coverage COVERAGE]
#                         [--deviation DEVIATION] [--threads THREADS]
# 
# Read Simulator Automation Wrapper for ART
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   --art ART             Path to ART executable
#   --samtools SAMTOOLS   Path to SAM tools executable and create BAM
#   --gzip GZIP           Path to gzip and compress all output
#   --fastadir FASTADIR   Process all FASTA files in directory
#   --fasta FASTA         Process individual FASTA files (comma delimited)
#   --output OUTPUT       Path to FASTQ/BAM output
#   --lengths LENGTHS     Read lengths (comma delimited)
#   --fragments FRAGMENTS Fragment length multipliers (comma delimited)
#   --coverage COVERAGE   Coverage amounts (comma delimited)
#   --deviation DEVIATION Standard deviation multipliers (comma delimited)
#   --threads THREADS     Number of threads to use
#
# Example: readSimulator.py --art C:\art\art_illumina.exe --fastadir C:\fastas --output C:\output --gzip C:\gzip\gzip.exe
#
# DEV TEST ARGUMENTS: --art Z:\Projects\Code\unh\mcbs913\art_bin_VanillaIceCream\art_illumina.exe --samtools Z:\Projects\Code\unh\mcbs913\samtools\samtools.exe --fastadir Z:\Projects\Code\unh\mcbs913\genomes --output Z:\Projects\Code\unh\mcbs913\temp --threads 1
 
import os
import subprocess
import argparse
import re
from concurrent.futures.thread import ThreadPoolExecutor

#-------------------- Arguments --------------------
argParser = argparse.ArgumentParser( description="Read Simulator Automation Wrapper for ART")
argParser.add_argument('--art', default="art_illumina.exe", help="Path to ART executable")
argParser.add_argument('--samtools', default="", help="Path to SAM tools executable and create BAM")
argParser.add_argument('--gzip', default="", help="Path to gzip and compress all output")
argParser.add_argument('--fastadir', default="", help="Process all FASTA files in directory")
argParser.add_argument('--fasta', default="", help="Process individual FASTA files (comma delimited)")
argParser.add_argument('--output', default='.', help="Path to FASTQ/BAM output")
argParser.add_argument('--lengths', default="50,150,250", help="Read lengths (comma delimited)")
argParser.add_argument('--fragments', default="2,3,4", help="Fragment length multipliers (comma delimited)")
argParser.add_argument('--coverage', default="10,100", help="Coverage amounts (comma delimited)")
argParser.add_argument('--deviation', default="0.1,0.4", help="Standard deviation multipliers (comma delimited)")
argParser.add_argument('--threads', default="4", help="Number of threads to use")

pairedArgument = ['', '-p']

# Build list of FASTA files to process
def buildFileList(passDirectory, passFiles):
    retFiles = []
    
    # If a directory was specified, add contained files to file list
    if passDirectory: 
        fastaRoot, fastaDirs, fastaFiles = next(os.walk(passDirectory))
        
        for currentFile in fastaFiles:
            retFiles.append({'full' : os.path.join(fastaRoot, currentFile), 'file' : currentFile}) 
        
    # If individual files were specified, add
    for currentFasta in passFiles.split(','):
        if currentFasta: retFiles.append({'full' : currentFasta, 'file' : os.path.basename(currentFasta)})
    
    return retFiles
    
# Output file prefix is constructed using the following convention: <taxon>-<se|pe>-<read length>-<fragment length>-<coverage>-[std deviation]  
def buildOutputPrefix(passOptions):
    # Obtain input filename without extension
    outputName = passOptions['outputDir'] + os.sep + os.path.splitext(passOptions['inputFile'])[0]

    # Append options
    outputName += "-{0}".format(['se', 'pe'][passOptions['paired']])    # SE/PE
    outputName += "-{0}".format(passOptions['optionRead'])              # Read length
    outputName += "-{0}".format(passOptions['optionFragment'])          # Fragment length
    outputName += "-{0}".format(passOptions['optionCoverage'])          # Coverage amount
    
    if passOptions['paired']:
        outputName += "-{0}".format(passOptions['optionDeviation'][3:]) # Standard deviation
        
    return outputName

# Sequence header is constructed from original, appended with nucleotide offset, and cleaned
def buildSequenceHeader(passLineFASTQ, passLineSAM, passPE, passStart):
    # Extract offset from SAM file 
    offsetMatch = re.match('^.*?\t.*?\t.*?\t(.*?)\t', passLineSAM)
    
    # Convert "|" delimiters to ":-:"        
    retHeader = passLineFASTQ.rstrip().replace('|', ':-:')
    
    # Temporarily remove pair ID 
    if passPE: retHeader = retHeader[:-2]
    
    # Append offset
    retHeader += ":-:{0}".format(offsetMatch.group(1))
    
    # Replace pair ID
    if passPE: retHeader += "/{0}".format(passStart + 1)

    retHeader += "\n"

    return retHeader

# In the case of paired-end reads, we clean up the filename by inserting a hyphen before the final number 
def cleanPairedName(passFilePrefix):
    os.rename("{0}1.fq".format(passFilePrefix), "{0}-1.fq".format(passFilePrefix))
    os.rename("{0}2.fq".format(passFilePrefix), "{0}-2.fq".format(passFilePrefix))

# Extract nucleotide offsets from SAM file and push into sequence description in FASTQ file
def transferOffsets(passPrefix, passPE, passStart):
    # Open both files for reading
    fastqName = "{0}.fq".format([passPrefix, "{0}-{1}".format(passPrefix, passStart + 1)][passPE])
    samName = "{0}.sam".format(passPrefix)
    
    fastqHandle = open(fastqName, 'r')
    samHandle = open(samName, 'r')
    
    # Open temporary file for writing (with .tmp extension)
    newName = "{0}.tmp".format(passPrefix)
    newHandle = open(newName, 'w')
    
    # Move to relevant part of SAM file
    for moveIdx in range(3 + passStart): samHandle.readline() 
    
    # Iterate through each line in FASTQ (without reading into memory)
    sectionIdx = 0
    for fastqLine in fastqHandle:
        # For header, generate and write new sequence header
        if sectionIdx == 0:
            header = buildSequenceHeader(fastqLine, samHandle.readline(), passPE, passStart)
            newHandle.write(header)
            
            # Skip subsequent SAM line for PE reads
            if (passPE): samHandle.readline()
        else:
            # Not a header, simply copy line
            newHandle.write(fastqLine)
            
        # Advance additional SAM entry for PE files, and FASTQ section index
        sectionIdx = (sectionIdx + 1) % 4

    # Close files
    fastqHandle.close()
    samHandle.close()
    newHandle.close()
    
    # Delete original FASTQ, replace with temporary file
    os.remove(fastqName)
    os.rename(newName, fastqName)

# Execute reads and associated functionality (ART, file cleanup, offset extraction, SAM->BAM)
def simulateReads(passOptions):
    # Build output filename and cache paired end option
    passOptions['outputPrefix'] = buildOutputPrefix(passOptions) 

    # Do not run Windows commands in a shell
    shellExec = (os.name != 'nt')
                        
    # Execute ART
    subprocess.call("{artExecutable} -i {inputFull} -l {optionRead} -m {optionFragment} -f {optionCoverage} {optionDeviation} -na -sam -o {outputPrefix}".format(**passOptions), shell=shellExec)
                        
    # Clean generated file names if paired end reads
    if passOptions['paired']: cleanPairedName(passOptions['outputPrefix'])
        
    # Transfer offsets from SAM into FASTQ
    if not passOptions['paired']: 
        transferOffsets(passOptions['outputPrefix'], False, 0)
    else:
        transferOffsets(passOptions['outputPrefix'], True, 0)
        transferOffsets(passOptions['outputPrefix'], True, 1)
        
    # Convert SAM to BAM
    if passOptions['samExecutable']: 
        subprocess.call("{0} view {1}.sam -S -b -o {1}.bam".format(passOptions['samExecutable'], passOptions['outputPrefix']), shell=shellExec)
    
    os.remove("{0}.sam".format(passOptions['outputPrefix']))
    
    # Compress files (FASTQ and BAM)
    if passOptions['gzipExecutable']:
        if not passOptions['paired']: 
            compressFiles = "{0}.fq ".format(passOptions['outputPrefix'])
        else:
            compressFiles = "{0}-1.fq {0}-2.fq ".format(passOptions['outputPrefix'])
        
        if passOptions['samExecutable']:
            compressFiles += "{0}.bam".format(passOptions['outputPrefix'])
            
        subprocess.call("{0} {1}".format(passOptions['gzipExecutable'], compressFiles), shell=shellExec)
    
#--------------------------- Main --------------------------------------
# Parse arguments
args = argParser.parse_args()

# Build the list of FASTA files to process
fastaFiles = buildFileList(args.fastadir, args.fasta)

# Execute for each file in FASTA Directory concurrently
with ThreadPoolExecutor(max_workers = int(args.threads)) as executor:
    for currentFile in fastaFiles:    
        # Nest for PE vs SE
        for optionPaired in pairedArgument:
            
            # Nest for read lengths
            for optionRead in args.lengths.split(','):
                
                # Nest for fragment lengths
                for optionFragment in args.fragments.split(','):
                    
                    # Nest for coverage amounts
                    for optionCoverage in args.coverage.split(','):
                        
                        # Nest for standard deviations
                        for optionDeviation in args.deviation.split(','):

                            # Calculate standard deviation
                            deviation = int(float(optionDeviation) * float(optionFragment) * float(optionRead))

                            # Build current options                    
                            currentOptions = {'artExecutable' : args.art,
                                              'samExecutable' : args.samtools,
                                              'gzipExecutable' : args.gzip,
                                              'outputDir' : args.output,
                                              'inputFile' : currentFile['file'],
                                              'inputFull' : currentFile['full'],
                                              'paired' : bool(optionPaired), 
                                              'optionPaired' : optionPaired, 
                                              'optionRead' : optionRead, 
                                              'optionFragment' : int(optionFragment) * int(optionRead),
                                              'optionCoverage' : optionCoverage,
                                              'optionDeviation' : ["-s {0}".format(deviation), ''][not optionPaired]
                                             }
    
                            # Simulate reads and all associated functionality with current options
                            executor.submit(simulateReads, currentOptions)
                            
                            # Only a single iteration for standard deviation is necessary for SE reads
                            if not optionPaired: break