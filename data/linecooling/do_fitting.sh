#! /bin/bash

mkdir tmp
for f in gamma_*.py
do
  python $f
done
