#!/usr/local/bin/python
#

import sys
import re
from optparse import OptionParser


def create_parser():
    """Creates the option parser used to parse through command-line args.
    
    Returns:
        options (like a dict, I guess): contains user-specified arguments.
        args (array): trailing arguments that are not paired with a flag.
    
    """
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)

    #SAM file argument
    parser.add_option("-s", "--sam", dest="sam",
                      help="read data from SAM file")

    #coords file argument
    parser.add_option("-c", "--coords", dest="coords",
                      help="read data from COORDS file")
    
    #mismatches argument (contained within the SNP file)
    parser.add_option("-m", "--mismatches", dest="mismatches",
                      help="read mismatch data from SNP file")
                      
                      
    options, args = parser.parse_args() 

    #raise error if required arguments are not supplied by the user
    if not options.mismatches:
        parser.error("missing --mismatches argument.")
    elif not options.coords:
        parser.error("missing --coords argument.")
    elif not options.sam:
        parser.error("missing --sam argument.")
    
    return options, args
   
def parse_coords_file(coords_file):
    """Parses through *.coords file to get contig names and alignments.
    
    Returns:
        contigs (dict): contains Contig objects representing contigs in the
                        assembly. Each one has alignment info to the reference.                        
    """
    contigs = dict()
    try:
        f = open(coords_file, 'r')
    except IOError as e:
        print >> sys.stderr, "coords file does not exist: " + coords_file
        exit(1)
        
    f.readline()
    f.readline()
    
    for line in f:
        elements = line.split("|")
        
        ref_start, ref_stop, = elements[0].rstrip().split(" ")
        name = elements[8].rstrip()
        name = re.sub(".*?(?=NODE)", '', name)
        c = Contig(name, ref_start, ref_stop)
        contigs[name] = c

    return contigs
    
def parse_sam_file(sam_file, contigs):
    """Parses through *.sam file to generate Read objects and add to contigs."""
    try:
        f = open(sam_file, 'r')
    except IOError as e:
        print >> sys.stderr, "SAM file does not exist: " + sam_file
        exit(1)
        
    for line in f:
        if line[0] == '@':
            continue
        elements = line.split("\t")
        read_name = elements[0]
        
        contig = contigs.get(elements[2])
        if not contig:
            continue
        
        aln_start = elements[3]
        read = Read(read_name, contig.ref_start, contig.ref_stop, aln_start)
        contig.reads.append(read)
        
        
def parse_mismatch_file(mismatch_file, contigs):
    """Parses through mismatches/indels found in snp file. Adds to contigs."""
    try:
        f = open(mismatch_file, 'r')
    except IOError as e:
        print >> sys.stderr, "snp file does not exist: " + mismatch_file
        exit(1)
        
   
    for line in f:
        elements = line.split("\t")
        contig = contigs.get(elements[1])
        if not contig:
            continue
            
        contig.errors.append(int(elements[5]))
           
class Contig(object):
    """Object that represents a contig in the assembly.
    
    Attributes:
        name (str): name of the contig.
        ref_start (int): start of the contig's alignment to the reference.
        ref_stop (int): stop of the contig's alignment to the reference.
        reads (list): contains Read objects that align to this contig.
        errors (list): list of errors in this contig, based on snp calls.         
    """         
    def __init__(self, name, ref_start, ref_stop):
        self.name = name
        self.ref_start = int(ref_start)
        self.ref_stop = int(ref_stop)
        self.reads = []
        self.errors = []
        
        
    def check_for_bad_alignment(self, verbose=False):
        """Iterates over each Read in read list. 
        
            If the Read is not aligned properly, do something...
            
        """ 
        bad_reads = 0
        good_reads = 0
        for read in self.reads:
            
            #is the read's start ref near our contig's start ref?
            read_start = read.read_ref_start
            contig_start = self.ref_start
            near_start = abs( read_start - contig_start ) <= 1000
            
            #if not, maybe it's near the end ref? (i.e. reverse direction)
            contig_stop = self.ref_stop
            near_stop = abs( read_start - contig_stop ) <= 1000
            
            #finally, maybe the read is just contained in the contig boundaries
            inside_boundaries = (read_start >= contig_start and 
                                 read_start <= contig_stop)
            
            #if neither is true, we know the read aligned somewhere wrong.
            if not (near_stop or near_start or inside_boundaries):
                
                #this is where we report the problem. Insert stuff here
                bad_reads += 1;
                read.correct_alignment = False
                
            #else, this read is fine
            else:
                good_reads += 1;     
        
        #don't give the user a report if we tell it to shut up.
        if not verbose:
            return
        print(self.name + " has " + str(good_reads) + " good reads and " + 
                str(bad_reads) + " bads reads.")
                
                
    def check_for_reads_near_errors(self, verbose=False):
        reads_near_errors = 0
        
        for read in self.reads:
            for error in self.errors:
                read_start = read.read_align_start
                if abs(read_start - error) >= 250:
                    continue
                
                reads_near_errors += 1
                read.near_error = True
                
        #don't give the user a report if we tell it to shut up.
        if not verbose:
            return
        print(self.name + " has " + str(reads_near_errors) + 
                " reads near errors.")
                
    def check_for_really_bad_reads(self, verbose=False):
        really_bad_reads = 0
        
        for read in self.reads:
            if read.near_error and not read.correct_alignment:
                really_bad_reads += 1
        
        #don't give the user a report if we tell it to shut up.
        if not verbose:
            return
        print(self.name + " has " + str(really_bad_reads) + 
                " reads that are near errors and are not aligned correctly." )
                

        
class Read(object):
    """Object that represents a read in the assembly.
    
    Attributes:
        name (str): name of the read (contains gi and positional information)
        contig_ref_start (int): start of where the contig that this read is a
                                part of aligns to the reference.
                                
        contig_ref_stop (int):  stop of where the contig that this read is a
                                part of aligns to the reference.
                                
        read_ref_start (int):   start of where the read was generated from
                                the reference.
                                
        read_ref_stop (int):    stop of where the read was generated from
                                the reference.
                                
        read_align_start (int): start of where the read aligns to the contig
                                that it is contained within.
                                
        read_align_stop (int): stop of where the read aligns to the contig
                                that it is contained within.
                                
        gi (str): the gi of the organism that this read came from
        correct_alignment(bool): True if, when this read aligns to a contig,
                                 that the contig's boundaries are near where
                                 the read SHOULD be in the reference.
                                 Defaults to: True
                                 
        near_error (bool):      True if this read is located near snps (errors)
                                in a contig. This may mean read is the source
                                of the mismatch/indel in the contig.
                                Defaults to: False
                                
                                
        
        
    """       
    def __init__(self, name, contig_ref_start, contig_ref_stop, aln_start):
        self.name = name
        self.contig_ref_start = int(contig_ref_start)
        self.contig_ref_stop = int(contig_ref_stop)
        self.read_align_start = int(aln_start)
        self.gi = None
        self.parse_name()
        self.correct_alignment = True
        self.near_error = False
        
    def parse_name(self):
        """Parses the read's name to get positional information."""
        m = re.match("gi:-:(\d*).*?:-:-(\d*):-:(\d*)", self.name)
        self.gi = m.groups()[0]
        self.read_ref_start = int(m.groups()[2])
        

#--------------------------- main --------------------------------------
def main():
    options, args = create_parser()
    contigs = parse_coords_file(options.coords)
    parse_sam_file(options.sam, contigs)
    parse_mismatch_file(options.mismatches, contigs)
    

    #loop over contigs and do something
    for key in contigs:
        contig = contigs[key]
        contig.check_for_bad_alignment(verbose=True)
        contig.check_for_reads_near_errors(verbose=True)
        contig.check_for_really_bad_reads(verbose=True)
        print("\n")

    
#execute main if this script was run via command-line
if __name__ == "__main__":
    main()
