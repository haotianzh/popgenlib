import numpy as np
import os
from ..base import Sample


def compute_ld_matrix(haplotype):
    matrix = haplotype.matrix
    rows, cols = matrix.shape
    freq_vec_1 = np.sum(matrix, axis=0) / rows
    freq_vec_0 = 1 - freq_vec_1
    freq_vec_1 = freq_vec_1[:, np.newaxis]
    freq_vec_0 = freq_vec_0[:, np.newaxis]
    product = np.dot(matrix.T, matrix) / rows
    p1q1 = np.dot(freq_vec_1, freq_vec_1.T)
    p0q0 = np.dot(freq_vec_0, freq_vec_0.T)
    ld = product - p1q1
    ld = ld ** 2 / (p1q1 * p0q0)
    return ld


def cluster(matrix):
    mat = matrix.copy()
    i = 0
    rows = matrix.shape[0]
    while i < rows:
        j = rows
        while j > i:
            temp = matrix[i:j, i:j]
            if np.sum(temp) / ((j - i) ** 2) > 0.1:  # threshold is 0.1 (can change here)
                mat[i:j, i:j] = 1
                i = j - 1
                break
            j = j - 1
        i = i + 1
    return mat


def compute_rh(data):
    rows, cols = data.shape
    remain = set()
    rh = 0
    for i in range(rows):
        if str(data[i, :]) not in remain:
            remain.add(str(data[i, :]))
            rh += 1
    rh = rh - cols
    return rh
