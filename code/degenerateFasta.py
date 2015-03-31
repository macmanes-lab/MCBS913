#!/usr/bin/python3
#
# Create a degenerate FASTA from original FASTA and GTF - degenerateFasta.py
#
# Toni Westbrook (anthonyw@wildcats.unh.edu)             
#
# usage: degenerateFasta.py [-h] [--fasta FASTA] [--gtf GTF] [--output OUTPUT]
#
# Create a degenerate FASTA from original FASTA and GTF
#
# optional arguments:
#   -h, --help       show this help message and exit
#   --fasta FASTA    Reference sequence (FASTA format)
#   --gtf GTF        Reference annotation (GTF / GFFv2 format)
#   --output OUTPUT  Output reference
# 
# --fasta Z:\temp\testing\PseudomonasFluorescensF113.fasta --gtf Z:\temp\testing\PseudomonasFluorescensF113.gff --output Z:\temp\testing\test.fasta

import argparse
import re

#-------------------- Arguments --------------------
argParser = argparse.ArgumentParser( description="Create a degenerate FASTA from original FASTA and GTF")
argParser.add_argument('--fasta', help="Reference sequence (FASTA format)")
argParser.add_argument('--gtf', help="Reference annotation (GTF / GFFv2 format)")
argParser.add_argument('--output', help="Output reference")

degenHash = { 'TTT' : 'TTY', 'TTC' : 'TTY', 'TTA' : 'TTR', 'TTG' : 'TTR',
              'CTT' : 'CTN', 'CTC' : 'CTN', 'CTA' : 'CTN', 'CTG' : 'CTN',
              'ATT' : 'ATH', 'ATC' : 'ATH', 'ATA' : 'ATH', 'ATG' : 'ATG',
              'GTT' : 'GTN', 'GTC' : 'GTN', 'GTA' : 'GTN', 'GTG' : 'GTN',
          
              'TCT' : 'TCN', 'TCC' : 'TCN', 'TCA' : 'TCN', 'TCG' : 'TCN',
              'CCT' : 'CCN', 'CCC' : 'CCN', 'CCA' : 'CCN', 'CCG' : 'CCN',
              'ACT' : 'ACN', 'ACC' : 'ACN', 'ACA' : 'ACN', 'ACG' : 'ACN',
              'GCT' : 'GCN', 'GCC' : 'GCN', 'GCA' : 'GCN', 'GCG' : 'GCN',              

              'TAT' : 'TAY', 'TAC' : 'TAY', 'TAA' : 'TAR', 'TAG' : 'TAR',
              'CAT' : 'CAY', 'CAC' : 'CAY', 'CAA' : 'CAR', 'CAG' : 'CAR',
              'AAT' : 'AAY', 'AAC' : 'AAY', 'AAA' : 'AAR', 'AAG' : 'AAR',
              'GAT' : 'GAY', 'GAC' : 'GAY', 'GAA' : 'GAR', 'GAG' : 'GAR',                           

              'TGT' : 'TGY', 'TGC' : 'TGY', 'TGA' : 'TGA', 'TGG' : 'TGG',
              'CGT' : 'CGN', 'CGC' : 'CGN', 'CGA' : 'CGN', 'CGG' : 'CGN',
              'AGT' : 'AGY', 'AGC' : 'AGY', 'AGA' : 'AGR', 'AGG' : 'AGR',
              'GGT' : 'GGN', 'GGC' : 'GGN', 'GGA' : 'GGN', 'GGG' : 'GGN',                           
            }

iupacHash = { 'A' : 'A', 'C' : 'C', 'G' : 'G', 'T' : 'T',
              'R' : 'AG', 'Y' : 'CT', 'S' : 'CG', 'W' : 'AT',
              'K' : 'GT', 'M' : 'AC', 'B' : 'CGT', 'D' : 'AGT',
              'H' : 'ACT', 'V' : 'ACG', 'N' : 'ACGT',
            }

iupacHashRev = dict (zip(iupacHash.values(), iupacHash.keys()))

class FastaFramework:
    """FASTA Framework for I/O and Parsing"""
        
    # Parse all headers and corresponding sequences from a FASTA file
    def __init__(self, passFileName):
        self.valid = False
        self.contigs = []
        self.selected = 0
        self.fileName = passFileName
                
        try:
            inFile = open( self.fileName, 'r' )
        except:
            print("Unable to open {0}".format(self.fileName))
            return

        # Iterate through FASTA file
        currentContigs = []
        for currentLine in inFile:
            if currentLine.startswith('>'):
                # Header - close out previous sequence, start new 
                if currentContigs: 
                    self.contigs[-1]['seq'] = ''.join(currentContigs)
                    self.contigs[-1]['dmask'] = bytearray(len(self.contigs[-1]['seq']))
                    
                self.contigs.append({'header':'', 'seq':'', 'dmask':''})
                self.contigs[-1]['header'] = currentLine.rstrip()
            else:
                # Sequence data
                currentContigs.append(currentLine.rstrip())

        # Close out final sequence and file, mark active
        if currentContigs: 
            self.contigs[-1]['seq'] = ''.join(currentContigs)
            self.contigs[-1]['dmask'] = bytearray(len(self.contigs[-1]['seq']))            
        inFile.close()
        self.valid = True

    # Get number of sequences contigs
    def getCount(self):
        return len(self.contigs)

    # Get sequence ID/description header for selected contig
    def getHeader(self):
        return self.contigs[self.selected]['header'][1:]

    # Get sequence for currently selected contig
    def getSequence(self):
        return self.contigs[self.selected]['seq']
    
    # Given two IUPAC codes, refine to the most ambiguous code that satisfies both constraints
    def refineDBase(self, passNew, passExist):
        # Check for existing mask values - if none present, use new nucleotide
        if not passExist:
            return ord(passNew)
        
        # Expand possible nucleotides of each IUPAC code
        newPossible = set(iupacHash[passNew])
        existPossible = set(iupacHash[chr(passExist)])
        
        # Refine ambiguity by calculating intersection
        refinePossible = newPossible.intersection(existPossible)
            
        # Sort lexicographically, then match back to IUPAC code
        refineNuc = ''.join(sorted(refinePossible))             
        
        return ord(iupacHashRev[refineNuc])
    
    # Given a CDS start and stop index, and direction, calculate the degenerate codon for each read frame     
    def calculateDMask(self, passStart, passEnd, passStrand, passFrame):
        # Generate iteration indices by strand direction
        if passStrand == '+':
            seqRange, step = range(passStart + passFrame, passEnd, 3), 1
        else:
            seqRange, step = range(passEnd - 1 - passFrame, passStart - 1, -3), -1
                        
        # Iterate through each codon in the coding sequence                   
        for seqIdx in seqRange:
            # Calculate degenerate codon
            codon = self.contigs[self.selected]['seq'][seqIdx:seqIdx + 3 * step:step]
            if codon not in degenHash.keys(): continue
            dCodon = degenHash[codon]

            # Iterate through each nucleotide in the codon
            for baseIdx in range(0, 3):
                # Calculate refined nucleotide
                dNucleotide = self.refineDBase(dCodon[baseIdx], self.contigs[self.selected]['dmask'][seqIdx + (baseIdx * step)])

                # Update Mask
                self.contigs[self.selected]['dmask'][seqIdx + (baseIdx * step)] = dNucleotide
    
    # Take the current degenerate mask, and apply to the nucleotide sequence where a value exists           
    def applyDMask(self):
        for seqIdx in range(0, len(self.contigs[self.selected]['seq'])):
            if not self.contigs[self.selected]['dmask'][seqIdx]: 
                self.contigs[self.selected]['dmask'][seqIdx] = ord(self.contigs[self.selected]['seq'][seqIdx])

        self.contigs[self.selected]['seq'] = self.contigs[self.selected]['dmask'].decode("ascii")
        
        
# Parse GTF line and return dictionary with relevant information
def getRowData(passLine):
    gtfColumns = ('seqname', 'source', 'feature', 'first', 'last', 'score', 'strand', 'frame', 'attribute')
        
    # Parse line
    gtfMatch = re.match('^(.*?)\\t(.*?)\\t(.*?)\\t(.*?)\\t(.*?)\\t(.*?)\\t(.*?)\\t(.*?)\\t(.*?)$', passLine)
    
    # Check for comment
    if not gtfMatch: return

    # Load into dictionary
    gtfDictionary = dict(zip(gtfColumns, gtfMatch.groups()))
    
    # Fixup default frames
    if gtfDictionary['frame'] == '.': gtfDictionary['frame'] = 0 
    
    return gtfDictionary

# Write output for each contig
def writeOutput(passFASTA, passFile):
    outputHandle = open(passFile, 'wb')
    
    for contigIdx in range(0, passFASTA.getCount()):
        passFASTA.selected = contigIdx
         
        # Write header
        outputHandle.write(bytes(">{0}\n".format(passFASTA.getHeader()), "ascii"))
         
        # Write output sequence
        sequenceBytes = bytes(passFASTA.getSequence(), "ascii")
        for sequenceIdx in range(0, len(sequenceBytes), 70):
            outputHandle.write(sequenceBytes[sequenceIdx:sequenceIdx+70] + bytes('\n', "ascii"))
     
    outputHandle.close()

#--------------------------- Main --------------------------------------
# Parse arguments
args = argParser.parse_args()

# Load FASTA reference
fastaRef = FastaFramework(args.fasta)
if not fastaRef.valid: exit() 

# Open GTF for reading
gtfHandle = open(args.gtf, 'r')

# Iterate through all GTF entries
for currentLine in gtfHandle:
    # Retrieve data for current GTF row
    currentData = getRowData(currentLine)
    
    # Skip until encountering a CDS
    if not currentData: continue
    if currentData['feature'] != 'CDS': continue;

    # Calculate the degenerate mask for the current CDS
    fastaRef.calculateDMask(int(currentData['first'])-1, int(currentData['last']), currentData['strand'], int(currentData['frame']))

# Close GTF file
gtfHandle.close()

# Apply the degenerate mask
fastaRef.applyDMask()

# Write output
writeOutput(fastaRef, args.output)