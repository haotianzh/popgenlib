from ..base import Haplotype
import os
import numpy as np

def sliding_windows(hap, window_size, step_size=None):
    if step_size is None:
        step_size = window_size
    windows = []
    positions, matrix = hap.positions, hap.matrix
    for i in range(0, hap.nsites, step_size):
        if i+window_size > hap.nsites:
            break
        window_mat = matrix[:, i:i+window_size]
        window_pos = positions[i:i+window_size]
        window_hap = Haplotype(matrix=window_mat, positions=window_pos)
        windows.append(window_hap)
    return windows


def write_vcf(trees, name):
    out = open(name, 'w')
    trees.write_vcf(out, ploidy=2)
    out.close()


def write_fasta(trees, name):
    out = open(name, 'w')
    trees.write_fasta(out)
    out.close()


def write_ms(trees, name):
    out = open(name, 'w')
    genotypes = trees.genotype_matrix()
    positions = []
    for v in trees.variants():
        positions.append(v.position)
    positions = np.array(positions).astype(int)
    out.write('ms {size} 1 -t {mu} -r {r} {length}\n'.format(size=100, mu=1e-8, r=1e-8, length=1e5))
    out.write('24263 40612 14324\n\n')
    out.write('//\n')
    out.write('segsites: {0}\n'.format(len(positions)))
    out.write('positions: ')
    out.write(' '.join(['{}'.format(val) for val in positions]))
    out.write('\n')
    for row in genotypes.T:
        out.write(''.join([str(val) for val in row]))
        out.write('\n')
    out.close()
