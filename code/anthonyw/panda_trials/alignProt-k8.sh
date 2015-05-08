#! /bin/sh

# Generate index and map
for readDir in /thomas1/mcbs913/anthony/panda_trials/prot/*
do
	for fullFile in $readDir/*/*.fna
	do
	   fileName=$(basename "$fullFile")
	   fullFileBase="${fullFile%%.*}"
	   echo "Aligning reads $readDir/reads.fq against $fullFileBase.fna and dumping to $fullFileBase-k8-o180.sam"
           /thomas1/mcbs913/anthony/mcbs913/bwa-fork/panda180 mem $fullFileBase.fna $readDir/reads.fq -t 20 -k 8 > $fullFileBase-k8-o180.sam
           cut -f1,2,3,4,5,6,7,8,9,10,11 $fullFileBase-k8-o180.sam | samtools view -Su - | samtools flagstat - > $fullFileBase-k8-o180.samstat
	done
done

