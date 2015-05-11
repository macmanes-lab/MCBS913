#!/usr/bin/python3
"""
Author: Junhong Chen

"""
import os
import argparse
import subprocess
import shlex


argParser = argparse.ArgumentParser( description="Novoalign")
argParser.add_argument('-i/--indir', dest = "indir", help="Dirctory containing the reads")
argParser.add_argument('-f/--ref', dest = "ref", help="Dirctory containing the references")
argParser.add_argument('-n/--ndx', dest = "ndx", help="Dirctory containing the ndx files")
argParser.add_argument('-o/--outdir', dest = "outdir", help="Output Directories for SAM files")


args = argParser.parse_args()

readpath = args.indir
if not readpath.endswith("/"):
    readpath = readpath + "/"
    
output = args.outdir
if not output.endswith("/"):
    output = output + "/"

ndx = args.ndx
if not ndx.endswith("/"):
    ndx = ndx + "/"
    
ref = args.ref

refdir = os.listdir(ref)

if not os.path.exists(args.outdir):
    print args.outdir,"not exists, make one!"
    os.makedirs(args.outdir)    
        
reads = set()    
    
for dirName, subdirList, fileList in os.walk(readpath):    
    for fname in fileList:
        tmp = fname.split(".")[0]
        reads.add(tmp[0:len(tmp)-1])
        

read_ref = dict()

for rd in reads:
    for rr in refdir:
        if rd.startswith(rr):
            read_ref[readpath+rd] = ref+"/"+rr
                    

for rdname in read_ref:
    rd1 = rdname+"1.fq.gz"
    rd2 = rdname+"2.fq.gz"
    for dirName, subdirList, fileList in os.walk(read_ref[rdname]):
        for fname in fileList:
            if fname.endswith("fna"):
                gname = fname.split(".")[0]
                ndxfile = [ndx + gname + ".ndx",ndx + gname + "_fd.ndx"]
                sam = [output + gname + ".sam",output + gname + "_fd.sam"]
                for i in range(2):
                    print "processing: ",ndxfile[i]
                    cmd = "python novoalignMT.py -nthreads 16 -d %s -f %s %s -i 1000,100 -o SAM -fname %s" % (ndxfile[i],rd1,rd2,sam[i]) 
                    subprocess.call(shlex.split(cmd))
                    
