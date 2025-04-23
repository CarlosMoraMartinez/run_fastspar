#!/bin/bash

INDIR=/media/DATOS22T/cmora/CORALS/inputdata_fastspar
OUTDIR=FastsparFull2
if [ -d "$OUTDIR" ]; then
    rm -rf "$OUTDIR"
fi

python call_fastspar.py -o $OUTDIR \
    -m mock \
    -a $INDIR/remove_tanda2_otus.tsv \
    -s age_class2 \
    --cleanup F \
    --nrand 10000 \
    --iterations 50 \
    --exclusion_iterations 10 \
    --exclusion_threshold 0.1 \
    --seed 123 \
    --iterations_parallel 10 \
    --threads 64


OUTDIR=FastsparAgeGroup2
if [ -d "$OUTDIR" ]; then
    rm -rf "$OUTDIR"
fi

python call_fastspar.py -o $OUTDIR \
    -m $INDIR/remove_tanda2_metad.tsv \
    -a $INDIR/remove_tanda2_otus.tsv \
    -s age_class2 \
    --cleanup F \
    --nrand 10000 \
    --iterations 50 \
    --exclusion_iterations 10 \
    --exclusion_threshold 0.1 \
    --seed 123 \
    --iterations_parallel 10 \
    --threads 64
