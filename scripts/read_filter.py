#!/usr/bin/env python3
"""
This script takes two inputs: a FASTQ file of reads and a PAF file of alignments to a reference.
It outputs (to stdout) reads from the FASTQ which pass this filter: at least half of the read's
bases must align to the reference. This aims to exclude reads which are total junk or are cross-
barcode contamination.
"""

import collections
import gzip
import sys


def main():
    read_filename = sys.argv[1]
    paf_filename = sys.argv[2]

    alignments = load_alignments(paf_filename)
    print(f'{len(alignments)} alignments loaded', file=sys.stderr)

    good_reads = set()
    for read, alignments in alignments.items():
        coverage = get_read_coverage(alignments)
        if coverage >= 0.5:
            good_reads.add(read)
    print(f'{len(good_reads)} reads pass filter', file=sys.stderr)

    print('Processing reads...', file=sys.stderr)
    for name, header, sequence, qualities in iterate_fastq(read_filename):
        if name in good_reads:
            print(header)
            print(sequence)
            print('+')
            print(qualities)


def load_alignments(paf_filename):
    alignments = collections.defaultdict(list)
    with open(paf_filename, 'rt') as f:
        for line in f:
            a = Alignment(line)
            alignments[a.read_name].append(a)
    return alignments


def get_read_coverage(alignments):
    """
    This function takes in all alignments for a read and it returns the fraction of the read which
    is covered by the alignments.
    """
    read_length = alignments[0].read_length
    bases = [0] * read_length
    for a in alignments:
        for i in range(a.read_start, a.read_end):
            bases[i] = 1
    return sum(bases) / read_length


class Alignment(object):

    def __init__(self, paf_line):
        line_parts = paf_line.strip().split('\t')
        if len(line_parts) < 11:
            sys.exit('Error: alignment file does not seem to be in PAF format')

        self.read_name = line_parts[0]
        self.read_length = int(line_parts[1])
        self.read_start = int(line_parts[2])
        self.read_end = int(line_parts[3])
        self.strand = line_parts[4]

        self.ref_name = line_parts[5]
        self.ref_length = int(line_parts[6]) / 2
        self.ref_start = int(line_parts[7])
        self.ref_end = int(line_parts[8])


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    https://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)
    unknown_file = open(str(filename), 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('\nError: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('\nError: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


def iterate_fastq(filename):
    with get_open_func(filename)(filename, 'rt') as fastq:
        for line in fastq:
            line = line.strip()
            if len(line) == 0:
                continue
            if not line.startswith('@'):
                continue
            name = line[1:].split()[0]
            header = line
            sequence = next(fastq).strip()
            _ = next(fastq)
            qualities = next(fastq).strip()
            yield name, header, sequence, qualities


if __name__ == '__main__':
    main()
