#!/bin/bash
for filename in metabat2_bins/*.sh; do
  sed -i 's/fna/fa/g' $filename
  msub $filename
done
