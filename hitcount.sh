#!/usr/bin/env bash

for cent in 1 2 11
do
  for k in $(seq 16 20)
  do
    echo "chr*_${k}mer_chr${cent}HOR_hit_test.gff"
    grep "kmer_hit" chr*_${k}mer_chr${cent}HOR_hit_test.gff | wc -l
  done
done
