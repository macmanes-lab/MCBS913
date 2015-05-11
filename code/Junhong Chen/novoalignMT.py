#!/usr/bin/python
"""
Author: Junhong Chen

"""
"""
A wrapper script for novoalign
"""
from multiprocessing import Process
import sys
import threading
import subprocess
import shlex
import time
from Bio import SeqIO
import gzip
from time import time

#options to run novoalign
cmd = "/opt/novocraft/novoalign"
tmpdir = "data_tmp"

#global variables
#ThreadLock = threading.Lock()
Threads = []


class ProcessParallel(object):
    """
    To Process the  functions parallely

    """    
    def __init__(self, jobs):
        """
        """
        self.jobs = jobs
        self.processes = []

    def fork_processes(self):
        """
        Creates the process objects for given function deligates
        """
        for job in self.jobs:
            proc  = Process(target=job[0],args = (job[1][0],job[1][1],job[1][2],))
            self.processes.append(proc)

    def start_all(self):
        """
        Starts the functions process all together.
        """
        for proc in self.processes:
            proc.start()

    def join_all(self):
        """
        Waits untill all the functions executed.
        """
        for proc in self.processes:
            proc.join()


def ParallelStart(jobs):
    procs =  ProcessParallel(jobs)
    procs.fork_processes()
    procs.start_all()
    procs.join_all()    


def initNovoJob():
    
    def callNovo(tid,cmd,tmp):
        with open(tmp+"out"+"_thread_"+str(tid)+".sam","w") as log:
            pro = subprocess.Popen(cmd,stdout = log)  
            pro.wait()
            
    info = []
    for i in range(1,options.nthreads+1):
        reads = [tmpdir+"/"+rd.split(".")[0].split("/")[-1]+"_thread_"+str(i)+".fq" for rd in options.reads]
        info.append([callNovo,[i,options.newCMD(reads),tmpdir+"/"]])   
        
    return info


class NovoThread (threading.Thread):
    
    def __init__(self, threadID, cmdline):
        
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.cmdLine = cmdline
        
    def run(self):
        print "Running thread: ",self.threadID
        args = self.cmdLine
        with open(tmpdir+"/"+"out"+"_thread_"+str(self.threadID)+".sam","w") as log:
            pro = subprocess.Popen(args,stdout = log)
            pro.wait()

    
        

class Options:
    
    def __init__(self):
        
        self.data = dict()

    def __getattr__(self,attr):
        
        return self.data[attr]
        
    def parse(self):
    
        ind = sys.argv.index("-nthreads")
        self.data["nthreads"] = int(sys.argv[ind+1])
        
        sys.argv.remove(sys.argv[ind+1])
        sys.argv.remove("novoalignMT.py")
        sys.argv.remove("-nthreads")
        
        if "-fname" in sys.argv:
            ind = sys.argv.index("-fname")
            self.data["outfile"] = sys.argv[ind+1]
            sys.argv.remove(sys.argv[ind+1])
            sys.argv.remove("-fname")
        else:
            self.data["outfile"] = None
        
        self.data["argv"] = [cmd] + sys.argv 
        
        self.data["reads"] = []
        ind = self.data["argv"].index("-f")
        self.data["reads"].append(self.data["argv"][ind+1])
        self.data["end"] = "se"
        if not self.data["argv"][ind+2].startswith("-"):
            self.data["reads"].append(self.data["argv"][ind+2])
            self.data["end"] = "pe"
            
    def newCMD(self,reads):
        
        ind = self.data["argv"].index("-f")
        ret = list(self.data["argv"])
        for i in range(1,len(reads)+1):
            ret[ind+i] = reads[i-1]
        return ret
 

def initReadings(fname,nthread,output):
    
    def batch_iterator(iterator, batch_size) :
        """Returns lists of length batch_size.
     
        This can be used on any iterator, for example to batch up
        SeqRecord objects from Bio.SeqIO.parse(...), or to batch
        Alignment objects from Bio.AlignIO.parse(...), or simply
        lines from a file handle.
     
        This is a generator function, and it returns lists of the
        entries from the supplied iterator.  Each list will have
        batch_size entries, although the final list may be shorter.
        """
        entry = True #Make sure we loop once
        while entry :
            batch = []
            while len(batch) < batch_size :
                try :
                    entry = iterator.next()
                except StopIteration :
                    entry = None
                if entry is None :
                    #End of file
                    break
                batch.append(entry)
            if batch :
                yield batch    
                
    
    tmpstr = fname.split(".")
                
    if tmpstr[-1] == "gz":
        hand = gzip.open(fname)
    elif tmpstr[-1] == "fq" or tmpstr[-1] == "fastq":
        hand = open(fname)
    else:
        print "Format for reading fils are not corect: ",fname
        sys.exit(-1)
        
    
                
    record_iter = SeqIO.parse(hand,"fastq")
    
    n = sum(1 for line in hand) / 4
    #print n
                
    hand = hand.seek(0)
                
    for i, batch in enumerate(batch_iterator(record_iter, int(n*1.0/nthread)+1)) :
        filename = output+fname.split(".")[0].split("/")[-1]+"_thread_%i.fq" % (i+1)
        handle = open(filename, "w")
        count = SeqIO.write(batch, handle, "fastq")
        handle.close()
        #print "Wrote %i records to %s" % (count, filename)    


        
        
        
def initThreads(options):
    
    info = []
    for i in range(1,options.nthreads+1):
        reads = [tmpdir+"/"+rd.split(".")[0]+"_thread_"+str(i)+".fq" for rd in options.reads]
        info.append([i,options.newCMD(reads)])

    #print info

    for elem in info:
        thread = NovoThread(elem[0], elem[1])
        Threads.append(thread)
    


def mainLock():
    
    for t in Threads:
        t.start()
    
    for t in Threads:
        t.join()    



def initTEMP(fname):
    
    subprocess.call("mkdir "+fname,shell=True)
    


def finaTEMP(fname):
    subprocess.call("rm -R "+fname,shell=True)

def mergeOutput(options):
    
    mg_cmd = "cat " + tmpdir+"/"+"out"+"_thread_1.sam"
    #remove header
    for i in range(2,options.nthreads+1):
        fname = tmpdir+"/"+"out"+"_thread_"+str(i)+".sam"
        out = tmpdir+"/"+"out"+"_thread_"+str(i)+"_noheader.sam"
        mg_cmd = mg_cmd + " " + out
        cmd = "samtools view -S " + fname
        args = shlex.split(cmd)
        with open(out,"w") as log:
            pro = subprocess.Popen(args,stdout = log,stderr = open("error_mgs_sam.merge","w"))
            pro.wait()

    subprocess.call("rm error_mgs_sam.merge",shell = True)
    #print mg_cmd
    #merge
    if options.outfile is None:
        outname = "out.sam"
    else:
        outname = options.outfile
        
    args = shlex.split(mg_cmd)
    with open(outname,"w") as log:
        pro = subprocess.Popen(args,stdout = log)    
        pro.wait()
   


if __name__ == "__main__":
    
    
    main_start = time()
    
    options = Options()
    options.parse()
    
    if options.outfile is not None and "/" in options.outfile:
        string = options.outfile.split("/")
        tmp = "/".join(string[0:len(string)-1])
        tmpdir = tmp + "/" + tmpdir
        
    #print options.reads
    #print options.argv
    #print options.newCMD(["reads1","reads2"])
    #print options.argv
    #print options.end
    
    initTEMP(tmpdir)
    
    init_start = time()
       
    print "Splitting reading files..........."
    if len(options.reads) == 1:
        initReadings(options.reads[0],options.nthreads,tmpdir+"/")
    else:
        jobs = [[initReadings, [options.reads[0],options.nthreads,tmpdir+"/"]],[initReadings, [options.reads[1],options.nthreads,tmpdir+"/"]]]
        ParallelStart(jobs)
    
    init_end = time()    
    
    run_start = time()
    
    ParallelStart(initNovoJob())
    
    #initThreads(options)
    #mainLock() 
    
    run_end = time()
    
    merge_start = time()
    
    print "Merging output files.............."
    mergeOutput(options)
    
    merge_end = time()
    
    
    main_end = time()    
    
    finaTEMP(tmpdir)
    

    
    print "Exiting Main Process"
    print "Total time: ",(main_end - main_start)/60.
    print "File Spliting time: ",(init_end - init_start)/60.
    print "running time: ",(run_end - run_start)/60.
    print "File merging time: ",(merge_end - merge_start)/60.
