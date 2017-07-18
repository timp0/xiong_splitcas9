#!/bin/bash
##Align MiSeq produced bisulfite data, but using very stringent settings to align to plasmid pool with differences of only 1-2 bases
##Especially for the spacing plasmid analysis datasets


#Ok - take in the root dir from command line arg 
##If command line arg empty - or not a dir - die                                      

if [ ! -d "$1" ]; then
    echo 
    echo "$1: No arg or dir doesn't exist"
    exit 1
else
    rawdir=$1
fi

##fastq pattern
if [ -z "$2" ]; then
    echo $2
    echo "No fastq pattern selected"
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


##If no geome, use bacteria as default

if [ -z "$4" ]; then
    genome="ecoli"
else
    genome=$4
fi


btidx=/mithril/Data/NGS/Reference/targeted/${genome}


scratchloc=/scratch/tmp/$RANDOM
echo $scratchloc

rm -R ${scratchloc}

rawread1=`ls ${rawdir}/${samp}*R1*.fastq.gz`
rawread2=`ls ${rawdir}/${samp}*R2*.fastq.gz`

mkdir -p ${scratchloc}

##Quality and adaptor trimming
~/Code/trim_galore_zip/trim_galore -q 30 --paired ${rawread1} ${rawread2} \
    -o ${scratchloc} &>${outdir}/${samp}_trimlog

##Change so that it has to be a perfect match to count as an alignment
~/Code/bismark/bismark --bam --bowtie2 \
    -p 4 \
    --most_valid_alignments 50 \
    --score_min L,0,0 \
    --non_directional \
    ${btidx} \
    -1 ${scratchloc}/*val_1.fq.gz \
    -2 ${scratchloc}/*val_2.fq.gz \
    --o ${outdir} \
    &>${outdir}/${samp}_galorepair_bismarklog


##Process aligned data
~/Code/bismark/bismark_methylation_extractor -p --multicore 4 --gzip --comprehensive \
    --genome_folder ${btidx} \
    ${outdir}/${samp}*bismark_bt2_pe.bam -o ${outdir} --no_header \
    &>${outdir}/${samp}_extractlog

~/Code/bismark/bismark2bedGraph --CX --dir ${outdir}/ --counts \
    -o ${samp}.bedGraph ${outdir}/C*${samp}*.txt.gz  \
    &>${outdir}/${samp}_bedlog

~/Code/bismark/coverage2cytosine -CX --genome_folder ${btidx} -o ${outdir}/${samp}.cyto.txt ${outdir}/${samp}.bismark.cov.gz

~/Code/bismark/bismark2report --alignment_report ${outdir}/${samp}*_PE_report.txt 


