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

#INDIR=/media/DATOS/cmora/2025_alejandra_pacbio/input_fastspar
#OUTDIR="/media/DATOS/cmora/2025_alejandra_pacbio/output_fastspar/Fastspar$varname"

#for varname in Group Weight
#do
#  OUTDIR="/media/DATOS/cmora/2025_alejandra_pacbio/output_fastspar/Fastspar$varname"
#  echo "$varname: $OUTDIR"
#  if [ -d "$OUTDIR" ]; then
#      rm -rf "$OUTDIR"
#  fi
#  
#  python call_fastspar.py -o $OUTDIR \
#      -m $INDIR/metadata.tsv \
#      -a $INDIR/filt_3samples.tsv \
#      -s $varname \
#      --cleanup F \
#      --nrand 10000 \
#      --iterations 50 \
#      --exclusion_iterations 10 \
#      --exclusion_threshold 0.1 \
#      --seed 123 \
#      --iterations_parallel 10 \
#      --threads 64
#done

INDIR=/media/DATOS/cmora/2025_alejandra_pacbio/input_fastspar2
OUTDIR1="/media/DATOS/cmora/2025_alejandra_pacbio/output_fastspar2/SIBO"
OUTDIR2="/media/DATOS/cmora/2025_alejandra_pacbio/output_fastspar2/Control"

python call_fastspar.py -o $OUTDIR1 \
      -m $INDIR/metadata_SIBO.tsv \
      -a $INDIR/filt_3samples_sibo.tsv \
      --cleanup T \
      --nrand 1000 \
      --iterations 50 \
      --exclusion_iterations 10 \
      --exclusion_threshold 0.1 \
      --seed 123 \
      --iterations_parallel 10 \
      --threads 64

python call_fastspar.py -o $OUTDIR2 \
      -m $INDIR/metadata_Control.tsv \
      -a $INDIR/filt_3samples_control.tsv \
      --cleanup T \
      --nrand 1000 \
      --iterations 50 \
      --exclusion_iterations 10 \
      --exclusion_threshold 0.1 \
      --seed 123 \
      --iterations_parallel 10 \
      --threads 64