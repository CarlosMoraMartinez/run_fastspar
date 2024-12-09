#!bin/bash

#Based on https://github.com/scwatts/fastspar
fastspar --yes --otu_table test/splitAbn_OTUs_FILT_T3_MSUS.csv \
    --correlation test/median_correlation_splitAbn_OTUs_FILT_T3_MSUS.tsv \
    --covariance test/median_covariance_splitAbn_OTUs_FILT_T3_MSUS.tsv \
    --iterations 5 -x 10 -e 0.1 -t 12 


fastspar_bootstrap --otu_table test/splitAbn_OTUs_FILT_T3_MSUS.csv \
    --number 20 \
    --prefix test/bootstrap_counts_splitAbn_OTUs_FILT_T3_MSUS/splitAbn_OTUs_FILT_T3_MSUS \
    -s 123 


parallel fastspar --yes --otu_table {} --correlation test/bootstrap_correlation_splitAbn_OTUs_FILT_T3_MSUS/cor_{/} --covariance test/bootstrap_correlation_splitAbn_OTUs_FILT_T3_MSUS/cov_{/} -i 5 ::: test/bootstrap_counts_splitAbn_OTUs_FILT_T3_MSUS/*

fastspar_pvalues --otu_table test/splitAbn_OTUs_FILT_T3_MSUS.csv --correlation test/median_correlation_splitAbn_OTUs_FILT_T3_MSUS.tsv \
    --prefix test/bootstrap_correlation_splitAbn_OTUs_FILT_T3_MSUS/cor_splitAbn_OTUs_FILT_T3_MSUS \
    --permutations 20 -t 12 \
    --outfile test/pvalues_splitAbn_OTUs_FILT_T3_MSUS.tsv