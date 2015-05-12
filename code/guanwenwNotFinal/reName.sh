#！/bin/bash
# Guanwen Wang
# rename some of the genome files

path='heterologous_references'
for dir in `ls $path`
do
    echo $dir
    cd $path/$dir        
    for files in `ls *`
    do
        mv $files `echo $files | sed "s/dir/$dir/g"`
        #mv $files `echo $files | sed 's/-/-$dir_/g'`
        #mv $files `echo $files | tr '-' '-$dir_'`
        #echo $files
    done
    cd ../..
done
