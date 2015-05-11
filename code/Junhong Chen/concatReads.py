"""
Author: Junhong Chen

"""
from Bio import SeqIO
import gzip
import sys
import os



pe1 = []
pe2 = []
pname = []


for dirName, subdirList, fileList in os.walk(sys.argv[1]):    
    for fname in fileList:
        tmp = fname.split(".")[0]
        tmp = tmp[:len(tmp)-1]
        if tmp not in pname:
            pname.append(tmp)
            pe1.append(dirName+"/"+tmp+"1.fq.gz")
            pe2.append(dirName+"/"+tmp+"2.fq.gz")
            
            


def concat(name,file_list):

    with open(name, 'w') as w_file:
        for filen in file_list:
            print 'working with',filen
            with gzip.open(filen, 'rU') as o_file:
                seq_records = SeqIO.parse(o_file, 'fastq')
                SeqIO.write(seq_records, w_file, 'fastq')
                
                
                
                
#print pe1
#print pe2
concat(sys.argv[2]+"-pe1.fq", pe1)
concat(sys.argv[2]+"-pe2.fq", pe2)
