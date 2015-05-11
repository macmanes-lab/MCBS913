#!/usr/bin/python3
"""
Author: Junhong Chen

"""
import os
import argparse
import subprocess

novoindex = "/opt/novocraft/novoindex"


argParser = argparse.ArgumentParser( description="Create nix files from novoindex")
argParser.add_argument('-i/--indir', dest = "indir", help="Dirctory containing the reference files")
argParser.add_argument('-o/--outdir', dest = "outdir", help="Output Directories")


args = argParser.parse_args()
    
if not os.path.exists(args.outdir):
    print args.outdir,"not exists, make one!"
    os.makedirs(args.outdir)    
    
# Set the directory you want to start from
rootDir = args.indir
for dirName, subdirList, fileList in os.walk(rootDir):
    print('Found directory: %s' % dirName)
    for fname in fileList:
        if fname.endswith("fna"):
            gname = fname.split(".")[0]
            fasta = dirName+"/"+fname
            ndx = args.outdir+"/"+gname+".ndx"
            cmd = "%s -t 16 %s %s" % (novoindex,ndx,fasta)
            subprocess.call(cmd,shell=True)
            
        if fname.endswith("fd"):
            gname = fname.split(".")[0]
            fasta = dirName+"/"+fname
            ndx = args.outdir+"/"+gname+"_fd.ndx"
            cmd = "%s -t 16 %s %s" % (novoindex,ndx,fasta)  
            subprocess.call(cmd,shell=True)
            
