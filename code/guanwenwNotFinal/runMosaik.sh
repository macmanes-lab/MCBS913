#!/bin/bash

#Guanwen Wang
#batch run Mosaik
#including MosaikAligner, MosaikBuild, MosaikJump

petail='-pe-250-1000-100-100'
setail='-se-50-100-10'

pe1fmt='-1.fq'
pe2fmt='-2.fq'
sefmt='.fq'

genomedir='/thomas1/mcbs913/junhong/dmapping/genomes/'
buildir='/thomas1/mcbs913/junhong/mosaik/temp/'
refdir='/thomas1/mcbs913/shared_data1/datasets/references/'
guandir='/thomas1/mcbs913/guanwen/MOSAIK-master/bin/'
outdir='/thomas1/mcbs913/junhong/mosaik/stats/'

declare -a gename=("HalobacillusHalophilusDSM2266" "MicrococcusLuteusNCTC2665" "AcidovoraxAvenaeATCC19860" "EscherichiaColiStrK-12SubstrMG1655" "PseudomonasFluorescensF113" "StaphylococcusEpidermidisATCC12228")
arraylength=${#gename[@]}

cd /opt/MOSAIK/bin

let count=0
declare -a fastdArr[24]
declare -a fastaArr[24]

echo "fastd"

for (( j=0; j<${arraylength}; j++ ));
do
    cd $refdir
    cd ${gename[$j]}

    for i in $( ls 1*fastd 2*fastd 3*fastd 4*fastd ); do

        fastdArr[$count]=$(echo $i | sed -r 's/.fastd/_degenRef/g')
        cd $guandir
  
        echo $buildir${gename[$j]}_bestread.dat
        echo $outdir${fastdArr[$count]}_alignedRes.dat
        echo $buildir${fastdArr[$count]}.dat
        echo $buildir${fastdArr[$count]}_jump_hs15      

#          ./MosaikAligner -in $buildir${gename[$j]}_bestread.dat -out $outdir${fastdArr[$count]}_alignedRes.dat -ia $buildir${fastdArr[$count]}.dat -hs 15 -mm 4 -j $buildir${fastdArr[$count]}_jump_hs15 -p 10 -a
#nnse /opt/MOSAIK/src/networkFile/2.1.26.se.100.005.ann -annpe /opt/MOSAIK/src/networkFile/2.1.26.pe.100.0065.ann   

#        ./MosaikBuild -fr $refdir${gename[$j]}/$i -oa $buildir${fastdArr[$count]}.dat &
        ./MosaikJump -ia $buildir${fastdArr[$count]}.dat -out $buildir${fastdArr[$count]}_jump_hs15 -hs 15 -iupac&
         

        echo $i
        echo ${fastdArr[$count]}
        ((count++))

    done
done

echo "fasta"
count=0

for (( j=0; j<${arraylength}; j++ ));
do
    cd $refdir
    cd ${gename[$j]}

    for i in $( ls 1*fasta 2*fasta 3*fasta 4*fasta ); do

        fastaArr[$count]=$(echo $i | sed -r 's/.fasta/_Normal/g')
        cd $guandir
        echo $buildir${gename[$j]}_bestread.dat
        echo $outdir${fastaArr[$count]}_alignedRes.dat
        echo $buildir${fastaArr[$count]}.dat
#        echo $buildir${fastdArr[$count]}_jump_hs15 
