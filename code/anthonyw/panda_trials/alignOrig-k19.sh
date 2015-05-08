#! /bin/sh

# Generate index and map
for readDir in /thomas1/mcbs913/anthony/panda_trials/data/orig/*
do
	for fullFile in $readDir/*/*.fna
	do
	   fileName=$(basename "$fullFile")
	   fullFileBase="${fullFile%%.*}"
	   echo "Aligning reads $readDir/reads.fq against $fullFileBase.fna and dumping to $fullFileBase-k19.sam"
           /usr/bin/bwa mem $fullFileBase.fna $readDir/reads.fq -t 20 -k 19 > $fullFileBase-k19.sam
           /usr/bin/samtools view $fullFileBase-k19.sam -Sbo $fullFileBase-k19.bam
           /usr/bin/samtools flagstat $fullFileBase-k19.bam > $fullFileBase-k19.samstat
	done
done

