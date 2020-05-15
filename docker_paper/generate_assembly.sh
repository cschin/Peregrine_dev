#!/bin/bash
echo $PWD/phased_reads/H1.fa > reads.lst
echo $PWD/phased_reads/unphased.fa >> reads.lst

echo yes | pg_run.py asm $PWD/reads.lst 4 4 4 4 4 4 4 4 4 --shimmer-r 3 \
            --with-consensus \
            --best_n_ovlp 8 \
            --output $PWD/asm-HG002-MHC-H1/

echo $PWD/phased_reads/H2.fa > reads.lst
echo $PWD/phased_reads/unphased.fa >> reads.lst

echo yes | pg_run.py asm $PWD/reads.lst 4 4 4 4 4 4 4 4 4 --shimmer-r 3 \
            --with-consensus \
            --best_n_ovlp 8 \
            --output $PWD/asm-HG002-MHC-H2/
