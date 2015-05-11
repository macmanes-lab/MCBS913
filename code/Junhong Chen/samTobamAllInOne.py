"""
Author: Junhong Chen

"""

import os
import sys
import subprocess

rootDir = sys.argv[1]


for dirName, subdirList, fileList in os.walk(rootDir):
    print('Found directory: %s' % dirName)
    for fname in fileList:
        if fname.endswith("sam"):
            print "sam file: ",fname
            filename = dirName + "/"+fname
            
            bname = dirName + "/"+fname.split(".")[0]
            oname = sys.argv[2] + "/"+fname.split(".")[0]
            #sname = sys.argv[3] + "/"+fname.split(".")[0]
            
            cmd = "samtools view -bS %s | samtools sort - %s" % (filename,oname)
            subprocess.Popen(cmd,shell=True)
            # cmd1 = "samtools flagstat %s > %s" % (oname+".bam",sname+".samstat")
            # subprocess.call(cmd1,shell=True)
