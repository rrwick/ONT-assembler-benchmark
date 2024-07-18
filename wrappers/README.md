# Assembler wrappers

This directory contains simple wrapper scripts for each assembly method, creating consistent inputs/outputs.

Each script takes the following three inputs:
1. the FASTQ file of ONT reads to be assembled
2. an estimated genome size (may or may not be used depending on the tools)
3. the sample name

Each script produces:
1. `sample__method.fasta`: the final assembly in FASTA format
3. `sample__method/`: a directory of all files made by the assembly method

For example, if I assembled sample `ABC` with the `flye-nano-hq` method, the final assembly would be named `ABC__flye-nano-hq.fasta` and there would be a `ABC__flye_nano-hq` directory containing the assembly files.

Requirements are listed in each script. All wrappers are set up to use 16 threads for steps that allow multithreading.
