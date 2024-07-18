#!/usr/bin/env bash

# Requirements:
#   Flye: https://github.com/mikolmogorov/Flye
#   Seqtk: https://github.com/lh3/seqtk

# Get inputs
reads=$1        # input reads FASTQ
genome_size=$2  # estimated genome size
sample=$3       # sample name

# Check that reads exist
if [[ ! -f "$reads" ]] ; then echo $reads" does not exist"; exit 1; fi

# Set up variables
threads=16
script=$(basename "$0")
method="${script%.*}"
output="$sample"__"$method"

# Do the assembly!
flye --nano-hq "$reads" --threads "$threads" --out-dir "$output"

# Copy the result into an uppercase single-line-per-sequence FASTA
seqtk seq -U "$output"/assembly.fasta > "$output".fasta
