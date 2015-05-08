#! /bin/bash

# Generate index and map
for lcv in {5..19}
do

   fileName=$(basename "$fullFile")
   fullFileBase="${fullFile%%.*}"
  
   echo current-$lcv 
   /thomas1/mcbs913/anthony/mcbs913/bwa-fork/panda mem seeds/2-Staphylococcus_aureus_JH1_NC_009632.fna seeds/reads.fq -t 20 -k $lcv > seeds/seedtest.sam
   cut -f1,2,3,4,5,6,7,8,9,10,11 seeds/seedtest.sam | samtools view -Su - | samtools flagstat - > seeds/seedtest-$lcv.samstat
done

