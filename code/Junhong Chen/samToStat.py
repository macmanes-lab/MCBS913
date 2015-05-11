"""
Author: Junhong Chen

"""
import os
import sys
import subprocess

rootDir = sys.argv[1]

i = 0
for dirName, subdirList, fileList in os.walk(rootDir):
    print('Found directory: %s' % dirName)
    for fname in fileList:
        if fname.endswith("bam"):
            print i,"bam file: ",fname
            filename = dirName + "/"+fname
            
            bname = dirName + "/"+fname.split(".")[0]
            oname = sys.argv[2] + "/"+fname.split(".")[0]
            
            cmd1 = "samtools flagstat %s > %s" % (bname+".bam",oname+".samstat")
            subprocess.call(cmd1,shell=True)
            i = i + 1
