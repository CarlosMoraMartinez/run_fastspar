#!/bin/bash

OUTDIR=test3
rm -rf $OUTDIR

python call_fastspar.py -o $OUTDIR \
    -m /home/carmoma/Documents/EarlyCause/results/descriptive1/METADATA_FILT.tsv \
    -a /home/carmoma/Documents/EarlyCause/results/descriptive1/OTUs_FILT.tsv \
    -s Time,Group \
    --cleanup T \
    --nrand 10 \
    --iterations 5 \
    --exclusion_iterations 10 \
    --exclusion_threshold 0.1 \
    --seed 123 \
    --iterations_parallel 10 \
    --permutations 10 \
    --threads 12
