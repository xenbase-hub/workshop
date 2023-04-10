#!/bin/bash


# Install SPAdes through bioconda
# $ conda install -c bioconda

# Ref. https://anaconda.org/bioconda/spades
# Ref. http://cab.spbu.ru/software/spades/

NUM_THREADS=8
DIR_FQ="../fastq.tx2/"

# Run trim+flash first.

FQ_e="TKLab202205mgi_XENTRtx_ctrlA_10m.extendedFrags.fastq"
FQ_nc="TKLab202205mgi_XENTRtx_ctrlA_10m.notCombined.fastq"
  
OUT=$(basename $FQ_e)
OUT=${OUT/.extendedFrags.fastq/}".spades"

if [ ! -e $OUT ]; then
    echo "Make "$OUT
    rnaspades.py --merged $FQ_e --12 $FQ_nc -t $NUM_THREADS -o $OUT
else
    echo $OUT" exists. Skip."
fi
