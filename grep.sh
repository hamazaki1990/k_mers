#!/usr/bin/env bash

for cent in 1 2 11
do
  for k in $(seq 20 22)
  do
    grep -v "##" chr*_${k}mer_chr${cent}HOR_hit_test.gff > ${k}mer_chr${cent}HOR_hit_test.gff
  done
done
