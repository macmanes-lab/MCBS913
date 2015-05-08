#! /bin/sh

# Generate index and map
for readDir in /thomas1/mcbs913/anthony/panda_trials/orig/*
do
	for fullFile in $readDir/*/*.fna
	do
	   fileName=$(basename "$fullFile")
	   fullFileBase="${fullFile%%.*}"
	   echo "Aligning reads $readDir/reads.fq against $fullFileBase.fna and dumping to $fullFileBase-k33.sam"
           /usr/bin/bwa mem $fullFileBase.fna $readDir/reads.fq -t 20 -k 33 > $fullFileBase-k33.sam
           /usr/bin/samtools view $fullFileBase-k33.sam -Sbo $fullFileBase-k33.bam
           /usr/bin/samtools flagstat $fullFileBase-k33.bam > $fullFileBase-k33.samstat
	done
done

