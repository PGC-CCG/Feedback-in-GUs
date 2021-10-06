#!/bin/bash

gradient_dir="/home/gfemer/GUs/Networks/geneInteractions_14-05-21_GRADIENT"
basedir=$(basename $gradient_dir)
echo -e $basedir
mkdir "GUs_${basedir}"

#screen -d -m -S Ptools bash -c 'pathway-tools -api -lisp'
#sleep 5

for ntw in $(ls ${gradient_dir}/GI-140521*.txt);
do
    base_error=$(basename $ntw | cut -f1 -d".")
    base=$(basename $ntw .txt)
    if [ $base != $base_error ]
    then
    echo -e "Removing GUs_${basedir}/${base_error} dir"
    rm -r "GUs_${basedir}/${base_error}"
    fi
    outdir="GUs_${basedir}/${base}/"
    if [ ! -d "$outdir" ]
    then 
    rm -r $outdir   
    mkdir $outdir
    echo -e "GUs for network: $base"
    nohup bash bin/pipeline_permissive_connectivity.sh $ntw $outdir > ${outdir}${base}.log 2>&1 &
    fi
done