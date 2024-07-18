#!/usr/bin/env python3
"""
This script takes a Unicycler GFA graph as input and it outputs (to stdout) an estimate of the
genome size. It does this by accounting for repeat contigs which need to be added to the total
multiple times.

This functionality is quite similar to Bandage's 'estimated sequence length' with one key
difference: it first splits the graph into connected components before normalising the depth. This
is to avoid inflated numbers from high-copy-number small plasmids. For example, a 3 kbp plasmid
that is in its own connected component should only add 3 kbp to the total, regardless of its depth.
"""

import collections
import sys


def main():
    gfa_filename = sys.argv[1]
    segments, links = load_graph(gfa_filename)
    components = split_graph(segments, links)
    estimated_genome_size = 0
    for segments in components:
        segments = normalise_depths(segments)
        for length, depth in segments.values():
            estimated_genome_size += length * depth
    print(estimated_genome_size)


def load_graph(gfa_filename):
    segments = {}
    links = collections.defaultdict(set)
    with open(gfa_filename) as f:
        for line in f:
            parts = line.rstrip().split('\t')
            if parts[0] == 'S':
                name, length, depth = get_segment_details(parts)
                segments[name] = (length, depth)
            if parts[0] == 'L':
                a, b = parts[1], parts[3]
                links[a].add(b)
                links[b].add(a)
    return segments, links


def get_segment_details(parts):
    name = parts[1]
    length = len(parts[2])
    depth = None
    for p in parts:
        if p.startswith('dp:f:'):
            depth = float(p[5:])
    if depth is None:
        sys.exit(f'Error: no depth found for segment {name}')
    return name, length, depth


def split_graph(segments, links):
    """
    Splits the graph into its connected components.
    """
    visited = set()
    components = []
    def dfs(node, current_component):
        stack = [node]
        while stack:
            v = stack.pop()
            if v not in visited:
                visited.add(v)
                current_component.add(v)
                stack.extend(links[v] - visited)
    for segment in segments:
        if segment not in visited:
            current_component = set()
            dfs(segment, current_component)
            component_segments = {seg: segments[seg] for seg in current_component}
            components.append(component_segments)
    return components


def normalise_depths(segments):
    depths_and_lengths = sorted((s[1], s[0]) for s in segments.values())
    median_depth = weighted_median(depths_and_lengths)
    new_segments = {}
    for name, length_and_depth in segments.items():
        length, depth = length_and_depth
        normalised_depth = int(round(depth / median_depth))
        new_segments[name] = (length, normalised_depth)
    return new_segments


def weighted_median(values_and_weights):
    total_weight = sum(x[1] for x in values_and_weights)
    target_weight = total_weight / 2
    weight_sum = 0
    for value, weight in values_and_weights:
        weight_sum += weight
        if weight_sum >= target_weight:
            assert value > 0.0
            return value
    assert False


if __name__ == '__main__':
    main()
