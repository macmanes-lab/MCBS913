#! /bin/sh

# Generate index and map
for fullFile in /thomas1/mcbs913/anthony/panda_trials/protnew/*/*/*.fna
do
   fileName=$(basename "$fullFile")
   fullFileBase="${fullFile%%.*}"
   /thomas1/mcbs913/anthony/mcbs913/bwa-fork/pandanew index $fullFileBase.fna $fullFileBase.gff
done

