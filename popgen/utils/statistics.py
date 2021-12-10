import numpy as np
import os
from ..base import Sample


def pairwise_ld(haplotype):
    """ Compute pairwise r^2 of LD.
    Input: a Haplotype object
    Return: a 2d numpy array
    """
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
    with np.errstate(divide='ignore', invalid='ignore'):
        ld = np.nan_to_num(ld ** 2 / (p1q1 * p0q0))
    return ld


def cluster_ld(matrix):
    """ Clustering based on precomputed pairwise LD matrix.
    Input: pairwise LD matrix (numpy.array)
    Return: clustering matrix (numpy.array)
    """
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


def rh(data):
    rows, cols = data.shape
    remain = set()
    rh = 0
    for i in range(rows):
        if str(data[i, :]) not in remain:
            remain.add(str(data[i, :]))
            rh += 1
    rh = rh - cols
    return rh


def pairwise_rf_distance():
    """ Compute pairwise Robinson-Foulds distance. """
    pass
