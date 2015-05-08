#!/usr/bin/python3
#
# Assembly automation via SPAdes - assemble.py
#
# Anthony Westbrook             
#
# usage: assemble.py [-h] [--spades SPADES] [--input INPUT] [--output OUTPUT]
#                    [--taxa TAXA] [--lengths LENGTHS] [--fragments FRAGMENTS]
#                    [--coverage COVERAGE] [--deviation DEVIATION]
#                    [--threads THREADS]
# 
# Assembly automation via SPAdes
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   --spades SPADES       Path to SPAdes python script
#   --input INPUT         Path to FASTQ files
#   --output OUTPUT       Path to assemblies
#   --taxa TAXA           Taxa (comma delimited)
#   --lengths LENGTHS     Read lengths (comma delimited)
#   --fragments FRAGMENTS Fragment length multipliers (comma delimited)
#   --coverage COVERAGE   Coverage amounts (comma delimited)
#   --deviation DEVIATION Standard deviation multipliers (comma delimited)
#   --threads THREADS     Number of threads to use
#
# Example: assemble.py --spades C:\spades\bin\spades.py --input C:\fastqs --output C:\output --taxa AcidovoraxAvenaeATCC19860
 
import os
import subprocess
import argparse
from concurrent.futures.thread import ThreadPoolExecutor

#-------------------- Arguments --------------------
argParser = argparse.ArgumentParser( description="Assembly automation via SPAdes")
argParser.add_argument('--spades', default="spades.py", help="Path to SPAdes python script")
argParser.add_argument('--input', default="", help="Path to FASTQ files")
argParser.add_argument('--output', default='.', help="Path to assemblies")
argParser.add_argument('--taxa', default="", help="Taxa (comma delimited)")
argParser.add_argument('--lengths', default="50,150,250", help="Read lengths (comma delimited)")
argParser.add_argument('--fragments', default="2,3,4", help="Fragment length multipliers (comma delimited)")
argParser.add_argument('--coverage', default="10,100", help="Coverage amounts (comma delimited)")
argParser.add_argument('--deviation', default="0.1,0.4", help="Standard deviation multipliers (comma delimited)")
argParser.add_argument('--threads', default="2", help="Number of threads to use")

pairedArgument = ['', '-p']

    
# Construct input prefix using the following convention: <taxa>-<se|pe>-<read length>-<fragment length>-<coverage>-[std deviation]  
def buildInputPrefix(passOptions):
    # Append options
    outputName = passOptions['taxa']                                    # Taxa
    outputName += "-{0}".format(['se', 'pe'][passOptions['paired']])    # SE/PE
    outputName += "-{0}".format(passOptions['optionRead'])              # Read length
    outputName += "-{0}".format(passOptions['optionFragment'])          # Fragment length
    outputName += "-{0}".format(passOptions['optionCoverage'])          # Coverage amount
    
    if passOptions['paired']:
        outputName += "-{0}".format(passOptions['optionDeviation'][3:]) # Standard deviation
        
    return outputName

# Execute reads and associated functionality (ART, file cleanup, offset extraction, SAM->BAM)
def assemble(passOptions):
    # Build input prefix and output directory
    inputPrefix = buildInputPrefix(passOptions)
    
    passOptions['inputFull'] = passOptions['inputDir'] + os.sep + inputPrefix 
    passOptions['outputFull'] = passOptions['outputDir'] + os.sep + inputPrefix
    
    # Do not run Windows commands in a shell
    shellExec = (os.name != 'nt')
                        
    # Execute SPAdes
    if not passOptions['paired']:
        #print("{spadesExecutable} -s1 {inputFull}.fq.gz -o {outputFull}".format(**passOptions))
        subprocess.call("{spadesExecutable} --s1 {inputFull}.fq.gz -o {outputFull}".format(**passOptions), shell=shellExec)
    else:
        #print("{spadesExecutable} -pe1-1 {inputFull}-1.fq.gz -pe1-2 {inputFull}-2.fq.gz -o {outputFull}".format(**passOptions))
        subprocess.call("{spadesExecutable} --pe1-1 {inputFull}-1.fq.gz --pe1-2 {inputFull}-2.fq.gz -o {outputFull}".format(**passOptions), shell=shellExec)

#--------------------------- Main --------------------------------------
# Parse arguments
args = argParser.parse_args()

# Assemble for each combination of prefix and argument
with ThreadPoolExecutor(max_workers = int(args.threads)) as executor:
    # Nest for taxa
    for optionTaxa in args.taxa.split(','):
        
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
                            currentOptions = {'spadesExecutable' : args.spades,
                                              'inputDir' : args.input,
                                              'outputDir' : args.output,
                                              'taxa' : optionTaxa,
                                              'paired' : bool(optionPaired), 
                                              'optionPaired' : optionPaired, 
                                              'optionRead' : optionRead, 
                                              'optionFragment' : int(optionFragment) * int(optionRead),
                                              'optionCoverage' : optionCoverage,
                                              'optionDeviation' : ["-s {0}".format(deviation), ''][not optionPaired]
                                             }
    
                            # Simulate reads and all associated functionality with current options
                            executor.submit(assemble, currentOptions)
                            
                            # Only a single iteration for standard deviation is necessary for SE reads
                            if not optionPaired: break