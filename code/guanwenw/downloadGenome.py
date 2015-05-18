#!/usr/bin/python3
import sys, getopt
import urllib.request
import os

'''
# Part of my team’s final submission
# Guanwen Wang
#
# run in terminal:
# ./downloadGenome.py -i list
# o == option
# a == argument passed to the o
#
# Note:list contains link like: you just need give any one of the links. 
# It will download .gff, .faa, .fna and store in the same directory.
#
# ftp://ftp.ncbi.nih.gov/genomes/Bacteria//Staphylococcus_pasteuri_SP1_uid226267/NC_022737.faa
# ftp://ftp.ncbi.nih.gov/genomes/Bacteria///Mycobacterium_smegmatis_JS623_uid184820/NC_019966.gff
#
'''

# Store input and output file names
ifile=''
ofile=''
inputFile=''

# Read command line args
myopts, args = getopt.getopt(sys.argv[1:],"i:o:")

for o, a in myopts:
    if o == '-i':
        ifile=a
    elif o == '-o':
        ofile=a
    else:
        print("Usage: %s -i input -o output" % sys.argv[0])

    # Display input and output file name passed as the args
    print ("Input file : %s and output file: %s" % (ifile,ofile) )

if ifile != "": 
    inputFile = ifile

i = 0 
lastDirName = ""
file_content = open(inputFile, 'r')
for line in file_content:
    token = line.split('\n')
    dirName = token[0].split('/')
    if dirName[6] == lastDirName:
        i = i + 1 
    else:
        lastDirName = dirName[6]
        i = 0 

    fileName = dirName[7].split('.')
    fileFaa = fileName[0]+".faa"
    fileFna = fileName[0]+".fna"
    fileGff = fileName[0]+".gff"

    ftp = token[0].split('.gff')
    urlPath = ftp[0]
    fullurl = [ftp[0]+".faa", ftp[0]+".fna", ftp[0]+".gff"]
    fileList = [fileFaa, fileFna, fileGff]

    dirList = [ dirName[6],'','']
    #print (fullurl)
    #print (fileList)

    if not os.path.exists( dirList[0] ):
        os.mkdir( dirList[0], 0o755 );

    for j in range(len(fileList)):
        path = dirList[0] + '/' + str(i + 1) + '-' + dirList[0] + '_' + fileList[j]

        print (path)
        print (fullurl[j])

        if not os.path.exists( path ):
            urllib.request.urlretrieve( fullurl[j], path)

        if os.path.exists( path ):
            print('Successfully downloaded\n')
        else:
            print('Failed Download\n')
