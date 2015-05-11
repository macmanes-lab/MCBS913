"""
Author: Junhong Chen

"""
import subprocess
import shlex
import sys

ndx = sys.argv[1]

rd1 = "/thomas1/mcbs913/junhong/metareads/metareads-pe1.fq"
rd2 = "/thomas1/mcbs913/junhong/metareads/metareads-pe2.fq"


cmd = "python novoalignMT.py -nthreads 16 -d %s -f %s %s -i 1000,100 -o SAM -fname %s" % (ndx+".ndx",rd1,rd2,"/thomas1/mcbs913/junhong/metareads/meta_standard.sam") 
subprocess.call(shlex.split(cmd))

cmd1 = "python novoalignMT.py -nthreads 16 -d %s -f %s %s -i 1000,100 -o SAM -fname %s" % (ndx+"_fd.ndx",rd1,rd2,"/thomas1/mcbs913/junhong/metareads/meta_degen.sam") 
subprocess.call(shlex.split(cmd1))
