#!/bin/bash

NUM_THREADS=4
MAX_OVERLAP=140


FQ1="../TKLab202205mgi_XENTRtx_ctrlA_10m_R1.raw.fastq.gz"
FQ2="../TKLab202205mgi_XENTRtx_ctrlA_10m_R2.raw.fastq.gz"

FQ2=${FQ1/_R1/_R2}
OUT_NAME=${FQ1/_R1.raw.fastq.gz/}
OUT_NAME=$(basename $OUT_NAME)

echo $FQ1 $FQ2 $OUT

flash $FQ1 $FQ2 -o $OUT_NAME -t $NUM_THREADS -M $MAX_OVERLAP --interleaved-output
