Anthony Westbrook <anthonyw@wildcats.unh.edu>
MCBS913 - Spring 2015

The following scripts and applications were developed for MCBS913.
Included are a sample micro.fasta and micro.gff (Micrococcus luteus) for
use with any program that requires a sequence and/or annotation file.

--------------------------------------------------------------------------------------------

Name: Protein Aligner (PALADIN)
Language: C (Make/gcc environment)
Description: This is the primary submission for the CS portion of the final project.

 (From main.c)

 PALADIN (Protein Alignment and Detection Interface)

 PALADIN is a protein sequence alignment tool based on the BWA source.  Like BWA, it aligns
 sequences via read-mapping using BWT.  PALADIN, however, offers the novel approach
 of aligning in the protein space.  During the index phase, it processes the reference genome's
 nucleotide sequences and GTF/GFF annotation containing CDS entries, first
 converting these transcripts into the corresponding protein sequences, then creating the BWT
 and suffix array from these proteins.  During the alignment phase, it attempts to find ORFs in 
 the read sequences, then converts these to protein sequences, and aligns to the reference 
 protein sequences. 

 PALADIN currently only supports single-end reads, and BWA-MEM based alignment.  It makes 
 use of many BWA parameters and is therefore compatible with many of its command line arguments.

 Filename: bwa-fork (directory - mixture of BWA source and custom code)

--------------------------------------------------------------------------------------------

Name: degenerateFasta
Language: Python3
Description: This is the primary submission for the CS portion of the final project.
 Create a degenerate FASTA from original FASTA and GTF.  Used for NovoAlign testing of degenerate sequence alignment.
 Error tested by Junhong and Guanwen using their scripts.
 
 usage: degenerateFasta.py [-h] [--fasta FASTA] [--gtf GTF] [--output OUTPUT]

 Create a degenerate FASTA from original FASTA and GTF

 optional arguments:
   -h, --help       show this help message and exit
   --fasta FASTA    Reference sequence (FASTA format)
   --gtf GTF        Reference annotation (GTF / GFFv2 format)
   --output OUTPUT  Output reference
 
 --fasta Z:\temp\testing\PseudomonasFluorescensF113.fasta --gtf Z:\temp\testing\PseudomonasFluorescensF113.gff --output Z:\temp\testing\degenerate.fasta
Filename: scripts/degenerateFasta.py


--------------------------------------------------------------------------------------------

Name: assemble
Language: Python3
Description: Assembly automation via SPAdes - runs a series of assemblies from ART generated reads 
 (arguments used in assemble.py should match those used with readSimulator.py)

 usage: assemble.py [-h] [--spades SPADES] [--input INPUT] [--output OUTPUT]
                    [--taxa TAXA] [--lengths LENGTHS] [--fragments FRAGMENTS]
                    [--coverage COVERAGE] [--deviation DEVIATION]
                    [--threads THREADS]
 
 Assembly automation via SPAdes
 
 optional arguments:
   -h, --help            show this help message and exit
   --spades SPADES       Path to SPAdes python script
   --input INPUT         Path to FASTQ files
   --output OUTPUT       Path to assemblies
   --taxa TAXA           Taxa (comma delimited)
   --lengths LENGTHS     Read lengths (comma delimited)
   --fragments FRAGMENTS Fragment length multipliers (comma delimited)
   --coverage COVERAGE   Coverage amounts (comma delimited)
   --deviation DEVIATION Standard deviation multipliers (comma delimited)
   --threads THREADS     Number of threads to use

 Example: assemble.py --spades C:\spades\bin\spades.py --input C:\fastqs --output C:\output --taxa AcidovoraxAvenaeATCC19860
Filename: scripts/assemble.py

--------------------------------------------------------------------------------------------
 
Name: convertHeaders
Language: Python3
Description: Convert FASTQ paired headers to include both alignments
 Will run for all files with a "-1.fq" and "-2.fq" suffix in the current directory
FileName: scripts/convertHeaders.py

--------------------------------------------------------------------------------------------

Name: readSimulator
Language: Python3
Description: Read Simulator via ART (Systematic creation of reads across a user specified range of parameters) - Codeveloped by Anthony Westbrook and Seth Hager
 usage: readSimulator.py [-h] [--art ART] [--samtools SAMTOOLS] [--gzip GZIP]
                         [--fastadir FASTADIR] [--fasta FASTA]
                         [--output OUTPUT] [--lengths LENGTHS]
                         [--fragments FRAGMENTS] [--coverage COVERAGE]
                         [--deviation DEVIATION] [--threads THREADS]
 
 Read Simulator Automation Wrapper for ART
 
 optional arguments:
   -h, --help            show this help message and exit
   --art ART             Path to ART executable
   --samtools SAMTOOLS   Path to SAM tools executable and create BAM
   --gzip GZIP           Path to gzip and compress all output
   --fastadir FASTADIR   Process all FASTA files in directory
   --fasta FASTA         Process individual FASTA files (comma delimited)
   --output OUTPUT       Path to FASTQ/BAM output
   --lengths LENGTHS     Read lengths (comma delimited)
   --fragments FRAGMENTS Fragment length multipliers (comma delimited)
   --coverage COVERAGE   Coverage amounts (comma delimited)
   --deviation DEVIATION Standard deviation multipliers (comma delimited)
   --threads THREADS     Number of threads to use

 Example: readSimulator.py --art C:\art\art_illumina.exe --fastadir C:\fastas --output C:\output --gzip C:\gzip\gzip.exe

 DEV TEST ARGUMENTS: --art Z:\Projects\Code\unh\mcbs913\art_bin_VanillaIceCream\art_illumina.exe --samtools Z:\Projects\Code\unh\mcbs913\samtools\samtools.exe --fastadir Z:\Projects\Code\unh\mcbs913\genomes --output Z:\Projects\Code\unh\mcbs913\temp --threads 1
Filename: scripts/readSimulatory.py

--------------------------------------------------------------------------------------------

Name: sam2csv
Language: Python3
Description: Convert all samtools flagstat outputs the current directory to a single CSV file
Filename: scripts/sam2csv.py

--------------------------------------------------------------------------------------------
Name: BWA/PALADIN test/trial scripts
Language: Bash
Description: Simple iteration scripts for running indexing and alignment runs across a number of parameters for both BWA  and PALADIN (Formerly known as PANDA)
Filename: panda_trials (directory)
