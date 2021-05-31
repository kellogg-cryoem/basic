#!/bin/bash

echo here

python ~/src/Rosetta/main/source/scripts/python/public/EnumerativeSampling/RunRosettaES.py \
 -rs runES.sh \
 -x RosettaES.xml \
 -f tnsc.fasta \
 -p input.pdb \
 -d map.mrc \
 -l $1 \
 -c 12 \
-n loop_$1 \
