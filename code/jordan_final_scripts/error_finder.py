#!/usr/local/bin/python
"""
    Script: error_finder v3.14
    Author: Jordan Ramsdell
    Desc:
        This script is meant to parse through metagenomic assemblie and find
        instances in which reads from different species map to the same contig.

        In such cases, the user can modify how these reads are filtered.

        The script can write the eronneous reads to a fasta file, and can
        also write the contigs in which these reads are located to a fasta file.
"""

import sys
import re
import cPickle
from os import path
import linecache
from optparse import OptionParser, OptionGroup, IndentedHelpFormatter
import mmap

# declare global variables and useful types
output = bytearray()
output_contig = bytearray()
counter = {}
total_counter = 0
contigs = {}
""":type : dict of [str, Contig]"""

def create_parser():
    """Creates the option parser used to parse through command-line args.
    
    Returns:
        options (like a dict, I guess): contains user-specified arguments.
        args (array): trailing arguments that are not paired with a flag.
    
    """
    usage = "usage: %prog --sam SAM --coords COORDS --snps SNPS [options]"

    help_formatter = IndentedHelpFormatter()
    help_formatter.set_long_opt_delimiter(" ")
    help_formatter.indent_increment = 0
    parser = OptionParser(usage, version="3.14",
                          formatter=help_formatter)

    group_required = OptionGroup(parser, "Required parameters")
    group_flags = OptionGroup(parser, "Filtering options")
    group_debug_flags = OptionGroup(parser, "Debug options")
    group_output = OptionGroup(parser, "Output options")

    # SAM file argument
    group_required.add_option("--sam", dest="sam",
                              help="Filepath to SAM file.")

    # coords file argument
    group_required.add_option("--coords", dest="coords",
                              help="Filepath to COORDS file.")

    # mismatches argument (contained within the SNP file)
    group_required.add_option("--snps", dest="mismatches",
                              help="Filepath to SNPS file.",
                              metavar="SNPS")

    # parameter to output contigs where misaligned reads are located
    parser.add_option("--fasta", dest="contigs_filename",
                      help="Fasta file containing contigs to extract"
                           " if reads mapping to the wrong "
                           "organisms are found.",
                      metavar="FASTA")

    # optional output name for contigs
    group_output.add_option("--contigs", dest="output_contigs",
                      help="Output name for contigs extracted from FASTA. "
                           "Requires --fasta to be set.", metavar="CONTIGS")

    # flag to output the reads that came from different species
    group_output.add_option("--out",
                      dest="output_reads", metavar="OUT",
                      help="When set, writes misaligned reads to OUT")

    # flag to only consider reads near snps
    group_flags.add_option("--near_errors", dest="near_error",
                           action="store_true",
                           help="When set, only considers reads that are "
                                "near (within 100 bp) mismatches found by "
                                "aligning the contigs to "
                                "the reference genomes via nucmer.")

    # flag to only consider reads near snps
    group_flags.add_option("--near_ends", dest="at_ends",
                           action="store_true",
                           help="When set, reads that are "
                                "aligned to the ends of contigs (within 250 bp)"
                                " are considered for checking for errors.")

    # option to filter contigs based on number of errors
    group_flags.add_option("--error_min", dest="error_min",
                           help="Specifies the minimum number of error reads"
                           "that must align to a contig before it is considered"
                           "an error that is worth reporting (default is 2).",
                           default=2, type="int", metavar="INT")

    # flag to show locations where reads from different species are mapping
    group_debug_flags.add_option("--show_regions", dest="show_regions",
                                 action="store_true",
                                 help="Shows the regions in contigs which "
                                      "reads from different species"
                                      " are mapping")

    # flag to show the misaligned read count from each species to each contig
    group_debug_flags.add_option("--show_counts", dest="read_count",
                                 action="store_true",
                                 help="Shows the number of reads in each contig"
                                      " that come from the wrong species.")

    parser.add_option_group(group_required)
    parser.add_option_group(group_output)
    parser.add_option_group(group_flags)
    parser.add_option_group(group_debug_flags)

    options, args = parser.parse_args()

    # raise error if required arguments are not supplied by the user
    error_message = ""
    if not options.mismatches:
        error_message = "missing --snps argument."
    elif not options.coords:
        error_message = "missing --coords argument."
    elif not options.sam:
        error_message = "missing --sam argument."

    if not error_message == "":
        parser.print_help()
        print "\nERROR: ", error_message
        exit(0)
        # parser.error(error_message)

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

        # PROBLEM: |'s in reference name are messing up parsing
        # SOLUTION: duct tape and bubblegum

        last_block = '%'.join(elements[4:])
        split = last_block.split(" ")

        ref_name, contig_name = split[1], split[2]
        ref_name = str(re.match('.*?%(.*?)%', ref_name).group(1))
        contig_name = contig_name.rstrip()
        ref_name = re.sub(".*?(?=NODE)", '', ref_name)

        c = Contig(ref_name, contig_name, ref_start, ref_stop)
        contigs[contig_name] = c

    f.close()

    return contigs


def parse_sam_file(sam_file, contigs):
    """Parses through *.sam file to generate Read objects and add to contigs."""
    try:
        f = open(sam_file, 'r+b')
    except IOError as e:
        print >> sys.stderr, "SAM file does not exist: " + sam_file
        exit(1)

    # FOR WINDOWS: mm = mmap.mmap(f.fileno(), 0, "ACCESS_READ")
    mm = mmap.mmap(f.fileno(), 0, mmap.MAP_PRIVATE, mmap.PROT_READ)
    line = mm.readline()
    idx = 0
    while line:

        if line[0] == '@':
            line = mm.readline()
            continue
        idx += 1
        elements = line.split("\t")

        contig = contigs.get(elements[2])

        if contig:
            read_name = elements[0]
            aln_start = elements[3]
            sequence = elements[9]
            contig.read_params.append([read_name, aln_start, sequence])
        line = mm.readline()

    print "The total number of reads in the sam file are: " + str(idx)
    mm.close()
    f.close()


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

    f.close()


class Contig(object):
    """Object that represents a contig in the assembly.
    
    Attributes:
        name (str): name of the contig.
        ref_name (str): name of the reference that this contig aligns to.
        ref_start (int): start of the contig's alignment to the reference.
        ref_stop (int): stop of the contig's alignment to the reference.
        reads (list): contains Read objects that align to this contig.
        errors (list): list of errors in this contig, based on snp calls.         
    """

    def __init__(self, ref_name, contig_name, ref_start, ref_stop):
        self.name = contig_name
        self.ref_name = ref_name
        self.ref_start = int(ref_start)
        self.ref_stop = int(ref_stop)
        self.reads = []
        self.errors = []
        self.sam_index = []
        self.read_params = []
        self.bad_read_count = 0
        self.is_bad_contig = False
        self.ranges = []
        m = re.match(".*?length_(.*?)_", self.name)
        self.length = int(m.group(1))

    def get_reads(self):
        for i in self.read_params:
            read_name = i[0]
            aln_start = i[1]
            sequence = i[2]
            self.reads.append(Read(read_name, self.ref_start, self.ref_stop,
                                   aln_start, sequence))

    def flush_reads(self):
        self.reads[:] = []
        self.read_params[:] = []

    def check_for_bad_alignment(self, verbose=False):
        """Iterates over each Read in read list. 
        
            If the Read is not aligned properly, do something...
            
        """
        bad_reads = 0
        good_reads = 0
        for read in self.reads:

            # is the read's start ref near our contig's start ref?
            read_start = read.read_ref_start
            contig_start = self.ref_start
            near_start = abs(read_start - contig_start) <= 1000

            # if not, maybe it's near the end ref? (i.e. reverse direction)
            contig_stop = self.ref_stop
            near_stop = abs(read_start - contig_stop) <= 1000

            # finally, maybe the read is just contained in the contig boundaries
            inside_boundaries = (contig_start <= read_start <= contig_stop)

            # if neither is true, we know the read aligned somewhere wrong.
            if not (near_stop or near_start or inside_boundaries):

                # this is where we report the problem. Insert stuff here
                bad_reads += 1
                read.correct_alignment = False

            # else, this read is fine
            else:
                good_reads += 1

        # don't give the user a report if we tell it to shut up.
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

        # don't give the user a report if we tell it to shut up.
        if not verbose:
            return
        print(self.name + " has " + str(reads_near_errors) +
              " reads near errors.")

    def check_for_really_bad_reads(self, verbose=False):
        really_bad_reads = 0

        for read in self.reads:
            if read.near_error and not read.correct_alignment:
                really_bad_reads += 1

        # don't give the user a report if we tell it to shut up.
        if not verbose:
            return
        print(self.name + " has " + str(really_bad_reads) +
              " reads that are near errors and are not aligned correctly.")

    def check_for_proper_alignment(self, options):
        start_offset = self.ref_start  # this is where the contig aligns to ref
        global output, counter, total_counter, output_contig
        potential_output = bytearray()

        for read in self.reads:
            # to figure out where read is mapping back to reference, we
            # add the contig's alignment start to the start of where
            # the read is aligning to the contig.
            actual_location = start_offset + read.read_align_start

            # we also know where the read SHOULD be in the reference by
            # looking in the read's header. We compare these two values.

            start_difference = abs(read.read_ref_start - actual_location)
            if options.at_ends:
                at_ends = read.read_align_start < 250 \
                          or abs(self.length - read.read_align_start) < 250
            else:
                at_ends = True

            if str(read.gi) != str(self.ref_name) and at_ends:

                # create a range for where contig aligns
                myname = [read.read_align_start, read.read_align_start + 250]
                merged = False
                for r in self.ranges:
                    if (r[0] <= myname[0] <= r[1]) \
                            or (r[1] >= myname[1] >= r[0]):
                        r[0] = myname[0] if myname[0] < r[0] else r[0]
                        r[1] = myname[1] if myname[1] > r[1] else r[1]
                        merged = True
                        break

                near_error = True if not options.near_error else read.near_error

                if not merged and near_error:
                    self.ranges.append(myname)
                    string = ">" + read.name + "___" + self.ref_name + \
                             "___" + self.name + "\n"

                    string += read.sequence + "\n"
                    potential_output.extend(string)

                total_counter += 1
                self.bad_read_count += 1
                if not counter.get(self.ref_name):
                    counter[self.ref_name] = {}
                if not counter[self.ref_name].get(self.name):
                    counter[self.ref_name][self.name] = {}
                if not counter[self.ref_name][self.name].get(read.gi):
                    counter[self.ref_name][self.name][read.gi] = 0

                counter[self.ref_name][self.name][read.gi] += 1

        if self.bad_read_count > options.error_min:
            output.extend(potential_output)
            self.is_bad_contig = True
            if options.show_regions and self.ranges:
                print self.name, "\n", self.ranges, "\n\n"


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
                                
        sequence (str):         The nucleotide sequence this read represents.  

        average_quality (float): The phred-score quality of the alignment.
                                
                                
        
        
    """
    def __init__(self, name, contig_ref_start, contig_ref_stop, aln_start,
                 sequence):
        self.name = name
        self.contig_ref_start = int(contig_ref_start)
        self.read_ref_start = 0
        self.contig_ref_stop = int(contig_ref_stop)
        self.read_align_start = int(aln_start)
        self.gi = None
        self.parse_name()
        self.correct_alignment = True
        self.near_error = False
        self.sequence = sequence

    def parse_name(self):
        """Parses the read's name to get positional information."""
        elements = re.compile(":-:").split(self.name)
        self.gi = str(elements[1])

        # more duct tape and bubblegum!
        if len(elements) == 4:
            self.read_ref_start = int(elements[3])
        else:
            self.read_ref_start = int(elements[5])


def find_errors(options):
    global contigs
    for key in contigs:
        contig = contigs[key]
        contig.get_reads()
        contig.check_for_reads_near_errors()
        contig.check_for_proper_alignment(options)
        contig.flush_reads()


def extract_bad_contigs():
    global contigs
    bad_contig_list = ""
    for key in contigs:
        contig = contigs[key]
        if contig.is_bad_contig:
            bad_contig_list += contig.name + "\n"
    myfile = open("bad_contig_list.txt", "w")
    myfile.write(bad_contig_list)
    myfile.close()


def extract_contigs(options):
    global contigs
    if not options.contigs_filename:
        return
    filename = options.contigs_filename
    if options.output_contigs:
        outname = options.output_contigs
    else:
        outname = "bad_contig_list.fasta"
    print "Writing contigs containing erroneous reads to " + outname
    contig_out = bytearray()
    with open(filename, "r") as f:
        file = f.read()
        for key in contigs:
            contig = contigs[key]
            if not contig.is_bad_contig:
                continue
            m = re.search(contig.name + "([^>]*)", file)
            if m:
                contig_out.extend(">" + contig.name + m.group(1))
            else:
                print "bad"

        myfile = open(outname, "w")
        myfile.write(contig_out)
        myfile.close()


# --------------------------- main --------------------------------------
def main():
    global contigs, output, total_counter
    options, args = create_parser()
    """:type : OptionParser"""

    contigs = parse_coords_file(options.coords)
    print "Parsing SAM file."
    parse_sam_file(options.sam, contigs)
    parse_mismatch_file(options.mismatches, contigs)
    print "Finding read errors."
    find_errors(options)

    # calculate total number of reads from wrong species
    if options.read_count:
        print "Showing contigs."
    final_total = 0
    for x in counter:
        for y in counter[x]:
            total = 0
            for z in counter[x][y]:
                total += counter[x][y][z]
            if total >= options.error_min:
                final_total += total
                if options.read_count:
                    print y, " -> ", counter[x][y], "\n\n"

    print "Number of erroneous reads before filtering:", total_counter
    print "Number of erroneous reads after filtering: " + str(final_total)

    # write outputs
    extract_contigs(options)
    out_name = "out.fa"
    if options.output_reads:
        out_name = options.output_reads
    print "Writing erroneous reads to " + out_name

    f = open(out_name, 'w')
    f.write(output)
    f.close()


# execute main if this script was run via command-line
if __name__ == "__main__":
    main()
