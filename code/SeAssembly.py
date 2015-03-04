#!/usr/bin/env python
__author__ = 'Guanwen Wang'

import sys, os, time, subprocess

SINFILE = "HHListS.txt"
PrintAll = 0

with open(SINFILE, "r") as ins:
    for line in ins:
        fileName = line.split('\n')
        dirArr = fileName[0].split('/')
        if PrintAll == 1:
            print "Read Line: %s" % (fileName[0])
        dirName = dirArr[6].split('.')
        if PrintAll == 1:
            print dirName

        path = dirName[0]
        os.mkdir( path, 0700 );
        time.sleep( 1 )
        cmd = "/opt/SPAdes-3.5.0-Linux/bin/spades.py -s " + fileName[0] + " -t 16 -o " + path + "/."
        if PrintAll == 1:
            print cmd 
        p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        while True:
            out = p.stderr.read(1)
            if out == '' and p.poll() != None:
                break
            if out != '':
                sys.stdout.write(out)
                sys.stdout.flush()
