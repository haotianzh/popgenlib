from ..base import Replicate, Haplotype
import os
import numpy as np
import scipy
from scipy import stats



def pairwise_ld(haplotype):
    """ Compute pairwise r^2 of LD.
    Input: a Haplotype object or a numpy ndarray
    Return: a 2d numpy array
    """
    if isinstance(haplotype, Haplotype):
        matrix = haplotype.matrix
    elif isinstance(haplotype, np.ndarray):
        matrix = haplotype
    else:
        Exception('Input should be an instance of either Haplotype or numpy.ndarray')
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


# def rh(data):
#     rows, cols = data.shape
#     remain = set()
#     rh = 0
#     for i in range(rows):
#         if str(data[i, :]) not in remain:
#             remain.add(str(data[i, :]))
#             rh += 1
#     rh = rh - cols
#     return rh


def stat(rates, pos, sequence_length, ne=1e5, window_size=50, step_size=50, bin_width=1e4, ploidy=1):
    snpsites = len(pos)
    rates = rates.reshape(-1)
    centers = []
    lens = []
    bounds = []
    for i in range(len(rates)):
        # take central point in each interval
        centers.append(pos[int(i*window_size+ step_size/2)])
        bounds.append((pos[i*window_size], pos[min(i*window_size+step_size, len(pos)-1)]))
        if i*window_size + step_size >= snpsites:
            last = len(rates)-1
        else:
            last = i*window_size + step_size
        lens.append(pos[last] - pos[i*window_size])
    lens = np.array(lens)
    scaledY = rates / lens / 2 / ploidy / ne
    v, bin_edges, _ = scipy.stats.binned_statistic(centers, scaledY, bins=sequence_length//bin_width) # range=(0,chrLength)
    return scaledY, bounds, [bin_edges[:-1], v]