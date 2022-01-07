import numpy as np
from ..base import Haplotype

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


