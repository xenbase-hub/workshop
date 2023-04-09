#!/bin/bash

## Make Tag Directory
makeTagDirectory tags.ctrl -format sam TKLab202203n_XENTRchip_ctrl_chr10.xenTro10_chr10.bam
makeTagDirectory tags.xtRFX2 -format sam TKLab202203n_XENTRchip_xtRFX2_chr10.xenTro10_chr10.sam 

## Find Peaks
findPeaks tags.xtRFX2 -style factor -o auto -i tags.ctrl

## Find Motifs
PEAKS="./tags.xtRFX2/peaks.txt"
GENOME="../data/genome/xenTro10.chr10.fa"
MOTIF_DIR="./motifs.xtRFX2"

findMotifsGenome.pl $PEAKS $GENOME $MOTIF_DIR

cat $MOTIF_DIR/knownMotifs/known?.motif > known_motifs.txt

## Annotate motifs
findMotifsGenome.pl $PEAKS $GENOME $MOTIF_DIR -find known_motifs.txt > rfx2.motif.pos.txt

GTF="../data/genome/XENTR_ncbi104.chr10.XB2023_04.gtf"
annotatePeaks.pl $PEAKS $GENOME -gtf $GTF -m known_motifs.txt > rfx2.motif.annot.txt
