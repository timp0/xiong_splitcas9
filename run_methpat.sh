#!/bin/bash

##Run methpat analysis on bismark aligned data - process starting from bam files

if [ ! -d "$1" ]; then
    echo 
    echo "$1: No arg or dir doesn't exist"
    exit 1
else
    rawdir=$1
fi

##bam pattern
if [ -z "$2" ]; then
    echo $2
    echo "No bam selected"
    exit 1
else 
    samp=$2
fi

##results directory - if empty or not a dir - die
if [ ! -d "$3" ]; then
    echo 
    echo "$3: No arg or dir doesn't exist - making"
    mkdir -p $3
    outdir=$3
    #exit 1
else
    outdir=$3
fi

##regfile 
if [ -z "$4" ]; then
    echo $4
    echo "No bam selected"
    exit 1
else 
    regions=$4
fi



echo ${rawdir}

mkdir -p ${outdir}/${samp}
    
bammy=`ls ${rawdir}${samp}/${samp}*.bam`
    
echo $bammy
    
    
~/Code/bismark/bismark_methylation_extractor ${bammy} --multicore 4 -o ${outdir}/${samp}/

methpat --amplicons ${regions} ${outdir}/${samp}/CpG_OT_${samp}*.txt --webassets online --logfile ${outdir}/${samp}.top.methpat.log \
    --html ${outdir}/${samp}.top.methpat.html >${outdir}/${samp}.top.methpat.tsv 

~/Code/methpat/methpat/methpat.py --filterpartial --amplicons ${regions} ${outdir}/${samp}/CpG_OT_${samp}*.txt --webassets online --logfile ${outdir}/${samp}.top.pf.methpat.log \
    --html ${outdir}/${samp}.top.pf.methpat.html >${outdir}/${samp}.top.pf.methpat.tsv 

#--min_cpg_percent PERCENT

methpat --amplicons ${regions} ${outdir}/${samp}/CpG_OB_${samp}*.txt --webassets online --logfile ${outdir}/${samp}.bot.methpat.log \
    --html ${outdir}/${samp}.bot.methpat.html >${outdir}/${samp}.bot.methpat.tsv 

~/Code/methpat/methpat/methpat.py --filterpartial --amplicons ${regions} ${outdir}/${samp}/CpG_OB_${samp}*.txt --webassets online --logfile ${outdir}/${samp}.bot.pf.methpat.log \
    --min_cpg_percent 1 \
    --html ${outdir}//${samp}.bot.pf.methpat.html >${outdir}/${samp}.bot.pf.methpat.tsv 


	


