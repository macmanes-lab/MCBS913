#! /bin/sh

# Generate index and map
for fullFile in /thomas1/mcbs913/anthony/panda_trials/data/prot/*/*/*.fna
do
   fileName=$(basename "$fullFile")
   fullFileBase="${fullFile%%.*}"
   /thomas1/mcbs913/anthony/mcbs913/bwa-fork/panda index $fullFileBase.fna $fullFileBase.gff
done

