#! /bin/sh

# Generate index and map
for fullFile in /thomas1/mcbs913/anthony/panda_trials/data/orig/*/*/*.fna
do
   fileName=$(basename "$fullFile")
   fileBase="${fileName%%.*}"
   /usr/bin/bwa index $fullFile
done

