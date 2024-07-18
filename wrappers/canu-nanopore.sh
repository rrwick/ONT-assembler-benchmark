#!/usr/bin/env bash

# Requirements:
#   Canu: https://github.com/marbl/canu
#   Seqtk: https://github.com/lh3/seqtk
#   canu_trim.py: https://github.com/rrwick/Trycycler/blob/main/scripts/canu_trim.py

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
canu -p canu -d "$output" -fast genomeSize="$genome_size" useGrid=false maxThreads="$threads" -nanopore "$reads"
canu_trim.py "$output"/canu.contigs.fasta > "$output"/canu.contigs.trimmed.fasta

# Copy the result into an uppercase single-line-per-sequence FASTA
seqtk seq -U "$output"/canu.contigs.trimmed.fasta > "$output".fasta
