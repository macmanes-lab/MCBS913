#! /bin/sh

# Generate index and map
for readDir in /thomas1/mcbs913/anthony/panda_trials/protnew/*
do
	for fullFile in $readDir/*/*.fna
	do
	   fileName=$(basename "$fullFile")
	   fullFileBase="${fullFile%%.*}"
	   echo "Aligning reads $readDir/reads.fq against $fullFileBase.fna and dumping to $fullFileBase-k9-o180.sam"
           /thomas1/mcbs913/anthony/mcbs913/bwa-fork/pandanew mem $fullFileBase.fna $readDir/reads.fq -t 15 -k 9 > $fullFileBase-k9-o180.sam
           cut -f1,2,3,4,5,6,7,8,9,10,11 $fullFileBase-k9-o180.sam | samtools view -Su - | samtools flagstat - > $fullFileBase-k9-o180.samstat
	done
done

