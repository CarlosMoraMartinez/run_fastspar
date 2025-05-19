#!/bin/bash

INDIR=/media/DATOS22T/cmora/CORALS/inputdata_fastspar
#OUTDIR=/media/DATOS22T/cmora/CORALS/output_fastspar/FastsparFull2
#if [ -d "$OUTDIR" ]; then
#    rm -rf "$OUTDIR"
#fi
#
#python call_fastspar.py -o $OUTDIR \
#    -m mock \
#    -a $INDIR/remove_tanda2_otus.tsv \
#    -s age_class2 \
#    --cleanup F \
#    --nrand 10000 \
#    --iterations 50 \
#    --exclusion_iterations 10 \
#    --exclusion_threshold 0.1 \
#    --seed 123 \
#    --iterations_parallel 10 \
#    --threads 64
#
#
#OUTDIR=/media/DATOS22T/cmora/CORALS/output_fastspar/FastsparAgeGroup2
#if [ -d "$OUTDIR" ]; then
#    rm -rf "$OUTDIR"
#fi
#
#python call_fastspar.py -o $OUTDIR \
#    -m $INDIR/remove_tanda2_metad.tsv \
#    -a $INDIR/remove_tanda2_otus.tsv \
#    -s age_class2 \
#    --cleanup F \
#    --nrand 10000 \
#    --iterations 50 \
#    --exclusion_iterations 10 \
#    --exclusion_threshold 0.1 \
#    --seed 123 \
#    --iterations_parallel 10 \
#    --threads 64


for varname in Sex Category_T0 hospital educ_m_discrete af_extraesc_m_00_Cat status_c2
do
  OUTDIR="/media/DATOS22T/cmora/CORALS/output_fastspar/Fastspar$varname"
  echo "$varname: $OUTDIR"
  if [ -d "$OUTDIR" ]; then
      rm -rf "$OUTDIR"
  fi
  
  python call_fastspar.py -o $OUTDIR \
      -m $INDIR/remove_tanda2_metad2.tsv \
      -a $INDIR/remove_tanda2_otus.tsv \
      -s $varname \
      --cleanup F \
      --nrand 10000 \
      --iterations 50 \
      --exclusion_iterations 10 \
      --exclusion_threshold 0.1 \
      --seed 123 \
      --iterations_parallel 10 \
      --threads 64
done